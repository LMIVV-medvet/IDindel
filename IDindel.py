import re
import argparse
import sys
import os

#Identify indels in an alignment of two sequences and reports the impacted annotations based of a GFF file of one of the sequences.
#usage: python IDIndel.py [-h] --fasta <FASTA file> --gff3 <GFF3 file> --out <BED output prefix>

#############################################################################
def parse_fasta(fasta_file):
    """
    Parses a FASTA file and returns a dictionary of sequences.
    
    Args:
        fasta_file (str): Input FASTA file.
    
    Returns:
        dict: Dictionary with sequence names as keys and sequences as values.
    """
    with open(fasta_file, 'r') as file:
        sequences = {}
        current_seq = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                current_seq = line[1:]
                sequences[current_seq] = ""
            else:
                sequences[current_seq] += line
                
    if len(sequences) > 2:
        sys.exit("Error: There are more than 2 sequences aligned in the fasta file.")

    return sequences

#############################################################################
def find_indels(seq1, seq2):
    """
    Finds insertions and deletions between two sequences.
    
    Args:
        seq1 (str): Reference sequence.
        seq2 (str): Second sequence.
    
    Returns:
        tuple: Two lists containing the start and end positions of indels in seq1 and seq2.
    """
    indels1 = []
    indels2 = []
    i = 0
    while i < len(seq1):
        if seq1[i] == '-':
            start = i
            while i < len(seq1) and seq1[i] == '-':
                i += 1
            indels2.append((start, i))
        elif seq2[i] == '-':
            start = i
            while i < len(seq2) and seq2[i] == '-':
                i += 1
            indels1.append((start, i))
        else:
            i += 1
    return indels1, indels2

#############################################################################
def write_combined_bed(indels1, indels2, output_file, seq1_name, seq2_name):
    """
    Writes a combined and sorted BED file with indels from two sequences.
    
    Args:
        indels1 (list): List of indels in the reference sequence.
        indels2 (list): List of indels in the second sequence.
        output_file (str): Output BED file.
        seq1_name (str): Name of the reference sequence.
        seq2_name (str): Name of the second sequence.
    """
    combined_indels = []
    for start, end in indels1:
        combined_indels.append([seq1_name, start, end])
    for start, end in indels2:
        combined_indels.append([seq2_name, start, end])
    combined_indels.sort(key=lambda x: int(x[1]))
    with open(output_file, 'w') as f:
        for indel in combined_indels:
            f.write('\t'.join(map(str, indel)) + '\n')

#############################################################################
def dealign_bed_file(input_bed_file, reference_sequence, prefix):
    """
    Parses a BED file with the indel coordinates and adjusts them so they correspond to the dealigned reference sequence.
    
    Args:
        input_bed_file (str): Input BED file.
        reference_sequence (str): Name of the reference sequence.
    """
    stored_length = 0
    previous_stored_length = 0
    adjusted_lines = []
    with open(input_bed_file, 'r') as file:
        lines = file.readlines()
    adjusted_lines = lines.copy()
    for i in range(len(lines)):
        line = lines[i].strip().split('\t')
        chrom, start, end = line[0], int(line[1]), int(line[2])
        if chrom != reference_sequence:
            interval_length = end - start
            stored_length += interval_length
            adjusted_lines[i] = f"{reference_sequence}\t{start - previous_stored_length}\t{start - previous_stored_length + 1}\t{interval_length}nt-insertion\n"
            previous_stored_length += interval_length
        else:
            adjusted_lines[i] = f"{chrom}\t{start-stored_length}\t{end-stored_length}\tdeletion\n"
    with open( prefix + '_adjusted_bed_file.bed', 'w') as output_file:
        output_file.writelines(adjusted_lines)

#############################################################################
def add_annotations_to_bed(bed_file, gff3_file, output):
    """
    Adds annotations from a GFF3 file to a BED file.
    
    Args:
        bed_file (str): Input BED file.
        gff3_file (str): Input GFF3 file.
    """
    annotations = {}
    with open(gff3_file, 'r') as file:
        for line in file:
            if not line.startswith("#"):
                parts = line.strip().split('\t')
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                annotation = parts[8]
                if chrom not in annotations:
                    annotations[chrom] = []
                annotations[chrom].append((start, end, annotation))
    with open(bed_file, 'r') as file:
        bed_lines = file.readlines()
    annotated_lines = []
    for bed_line in bed_lines:
        parts = bed_line.strip().split('\t')
        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        overlapping_annotations = []
        if chrom in annotations:
            for ann_start, ann_end, annotation in annotations[chrom]:
                if ann_start <= end and ann_end >= start:
                    overlapping_annotations.append(annotation)
        if overlapping_annotations:
            annotated_line = bed_line.strip() + '\t' + ';'.join(overlapping_annotations) + '\n'
        else: 
            annotated_line = bed_line.strip() + '\t' + "intergenic" + '\n'
        annotated_lines.append(annotated_line)
    with open(output, 'w') as output_file:
        output_file.writelines(annotated_lines)

#############################################################################
def find_chromosome_names(gff3_file):
    """
    Finds and returns the sequence name in a GFF3 file.
    Exits with an error if there are multiple sequence names.
    
    Args:
        gff3_file (str): Path to the GFF3 file.
    
    Returns:
        str: Sequence name.
    """
    chromosome_names = set()
    
    with open(gff3_file, 'r') as file:
        for line in file:
            if not line.startswith("#"):
                parts = line.strip().split('\t')
                chrom = parts[0]
                chromosome_names.add(chrom)
    
    if len(chromosome_names) > 1:
        sys.exit("Error: There are more than 1 sequence name in the GFF3 file.")
    
    return chromosome_names.pop()

#############################################################################
def main():
    # Get argurments
    parser = argparse.ArgumentParser(description="Identify indels in an alignment of two sequences and reports the impacted annotations based of a GFF file of one of the sequences.")
    parser.add_argument('-f','--fasta', metavar='fasta', required=True, help="Path to the input FASTA alignment file.")
    parser.add_argument('-g','--gff3', metavar='gff3', required=True, help="Path to the input GFF3 annotation file.")
    parser.add_argument('-o','--out', metavar='string', required=True, help="Path to the BED output file with prefix.")

    args = parser.parse_args()
    
    fasta_file = args.fasta
    gff3_file = args.gff3
    output_file = args.out + "-Indel-annotations.bed"
    
    
    if not os.path.exists(fasta_file):
        sys.exit("Error: the input alignement file is not found.")

    if not os.path.exists(gff3_file):
        sys.exit("Error: the input GFF3 file is not found.")

    
    prefix = fasta_file.rsplit('.', 1)[0]
    
    #Identify the reference sequence and the query
    reference_sequence = find_chromosome_names(gff3_file)
    sequences = parse_fasta(fasta_file)
    seq_names = list(sequences.keys())
    seq1_name = reference_sequence
    
    if not reference_sequence in seq_names:
        sys.exit("Error: The reference sequence was not found in the alignment file. Ensure that the sequence name in the GFF3 is the same as one of the aligned sequences.")
    
    if seq_names[0] != reference_sequence:
        seq2_name = seq_names[0]
    else:
        seq2_name = seq_names[1]
    seq1 = sequences[seq1_name]
    seq2 = sequences[seq2_name]
    
    #Make the bed file ith the indel positions
    indels1, indels2 = find_indels(seq1, seq2)
    combined_bed_file = prefix + '_combined_sorted_indels.bed'
    write_combined_bed(indels1, indels2, combined_bed_file, seq1_name, seq2_name)
    
    #Adjust the bed file to the dealigned sequence and add the annotations
    dealign_bed_file(combined_bed_file, reference_sequence, prefix)
    add_annotations_to_bed(prefix + '_adjusted_bed_file.bed', gff3_file, output_file)
    
    #Clean up
    os.remove( prefix + '_combined_sorted_indels.bed')
    os.remove( prefix + '_adjusted_bed_file.bed')

    print("IDindel ran successffully")

if __name__ == "__main__":
    main()





