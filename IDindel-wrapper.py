
#!/usr/bin/env python3
"""
Wrapper to run IDindel.py over all FASTA files in a directory, passing:
  python <path_to_wrapper>/IDindel.py -f <fasta_file> -g <GFF_file> -o <output_dir>/<prefix>.bed

Requirements:
  - GFF file must exist and end with .gff or .gff3.
  - FASTA files must end with .fa or .fasta and contain exactly 2 sequences.
"""

import argparse
import sys
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


def error_and_exit(message: str, code: int = 1):
    print(f"[ERROR] {message}")
    sys.exit(code)


def compute_prefix(path: Path) -> str:
    """Return filename without ANY extensions."""
    p = path
    for _ in p.suffixes:
        p = p.with_suffix("")
    return p.name


def find_fasta_files(directory: Path):
    """Return all FASTA files (.fa or .fasta) in the directory."""
    return sorted([f for f in directory.iterdir() if f.is_file() and f.suffix.lower() in [".fa", ".fasta"]])


def validate_gff_file(gff_path: Path) -> bool:
    return gff_path.exists() and gff_path.is_file() and gff_path.suffix.lower() in [".gff", ".gff3"]


def validate_fasta_content(file_path: Path) -> bool:
    """Ensure FASTA file contains exactly 2 sequences."""
    count = 0
    try:
        with file_path.open() as f:
            for line in f:
                if line.startswith(">"):
                    count += 1
    except Exception:
        return False
    return count == 2


def run_idindel(python_exe: str, id_script: Path, fasta_path: Path, gff_path: Path, out_file: Path):
    """Run IDindel.py with required arguments."""
    cmd = [
        python_exe, str(id_script),
        "-f", str(fasta_path),
        "-g", str(gff_path),
        "-o", str(out_file)
    ]
    print(f"[INFO] Running: {' '.join(cmd)}")
    completed = subprocess.run(cmd)
    return fasta_path, completed.returncode


def main():
    parser = argparse.ArgumentParser(description="Run IDindel.py on all FASTA files in a directory with a GFF file.")
    parser.add_argument("-d","--dir", required=True, help="Directory containing FASTA alignment files (.fa/.fasta).")
    parser.add_argument("-g","--gff3", required=True, help="Path to the reference GFF file (.gff or .gff3).")
    parser.add_argument("-o","--out", required=True, help="Output directory (created if missing).")
    parser.add_argument("-t","--threads", type=int, default=1, help="Number of threads to use (default: 1).")
    args = parser.parse_args()

    # Locate IDindel.py in the same directory as this wrapper
    wrapper_dir = Path(__file__).resolve().parent
    id_script = wrapper_dir / "IDindel.py"
    if not id_script.exists():
        error_and_exit(f"IDindel.py not found in wrapper directory: {id_script}")

    # Validate input directory
    input_dir = Path(args.dir)
    if not input_dir.exists() or not input_dir.is_dir():
        error_and_exit(f"Input directory does not exist or is not a directory: {input_dir}")

    # Validate GFF file
    gff_path = Path(args.gff3)
    if not validate_gff_file(gff_path):
        error_and_exit(f"Invalid GFF file: {gff_path}")

    # Validate workers
    if args.threads < 1:
        error_and_exit("--threads must be >= 1.")

    # Prepare output directory
    out_dir = Path(args.out)
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        error_and_exit(f"Cannot create output directory '{out_dir}': {e}")

    # Find FASTA files
    fasta_files = find_fasta_files(input_dir)
    if not fasta_files:
        error_and_exit("No FASTA files (.fa/.fasta) found in the specified directory.", code=0)

    print(f"[INFO] Found {len(fasta_files)} FASTA file(s). Threads: {args.threads}. Output dir: {out_dir}")

    # Validate FASTA content
    invalid_files = [f for f in fasta_files if not validate_fasta_content(f)]
    if invalid_files:
        print("[ERROR] The following files do NOT contain exactly 2 sequences:")
        for f in invalid_files:
            print(f"  - {f}")
        sys.exit(1)

    # Prepare tasks
    tasks = [(f, out_dir / f"{compute_prefix(f)}") for f in fasta_files]
    python_exe = sys.executable
    failures = 0

    # Run in parallel
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_file = {
            executor.submit(run_idindel, python_exe, id_script, f, gff_path, out_file): f
            for (f, out_file) in tasks
        }
        for future in as_completed(future_to_file):
            file_done = future_to_file[future]
            try:
                file_path, rc = future.result()
            except Exception as exc:
                print(f"[ERROR] Exception for file {file_done}: {exc}")
                failures += 1
            else:
                if rc != 0:
                    print(f"[ERROR] Command failed (rc={rc}) for file: {file_path}")
                    failures += 1
                else:
                    print(f"[OK] Completed: {file_path}")

    if failures:
        print(f"[DONE] Completed with {failures} failure(s).")
        sys.exit(1)
    else:
        print("[DONE] All runs completed successfully.")
        sys.exit(0)


if __name__ == "__main__":
    main()
