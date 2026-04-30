import os
import pandas as pd
import argparse

#Creates a presence/absence matrix of the indels in multiple bed files produced by the IDIndel.py script.
#usage: python MatrixBuilder.py [-h] --input <input directory> --prefix <output prefix> [--start_pos <INTEGER> --end_pos <INTEGER>]


############################################################
def concatenate_bed_files(input_directory, exclusions):
    # List to store dataframes
    dataframes = []

    # Iterate over all files in the directory
    for filename in os.listdir(input_directory):
        if filename.endswith("-Indel-annotations.bed"):

            #process if file is not empty due to no annotations
            if os.path.getsize(input_directory+"/"+filename) > 0:

                # Read the bed file into a dataframe
                df = pd.read_csv(os.path.join(input_directory, filename), sep='\t', header=None)

                # Modify the first column to be the filename without the suffix
                df[0] = filename.replace("-Indel-annotations.bed", "")

                # Append the dataframe to the list
                dataframes.append(df)

            else :
                sequence = filename.replace("-Indel-annotations.bed", "")
                append_line(exclusions, sequence)

    # Concatenate all dataframes
    concatenated_df = pd.concat(dataframes)
    return concatenated_df

##################################################################
def build_presence_absence_matrix(concatenated_df, exclusions, start_pos=None, end_pos=None):
    # Create an empty dictionary to store presence/absence data
    presence_absence_dict = {}
    exclusion_dict={}
    strain_count={}

    # Iterate over each row in the concatenated dataframe
    for index, row in concatenated_df.iterrows():
        strain = row[0]
        start = row[1]
        end = row[2]
        indel_type = row[3]

        # Create a dictionary to count the number of lines for a given strain
        if strain not in strain_count:
            strain_count[strain] = 1
        else:
            strain_count[strain] += 1 	

        # Create a dictionary to count the number of excluded lines for a given strain
        if strain not in exclusion_dict:
            exclusion_dict[strain] = 0

        # Apply the start and end position filters if provided and count the number of omitted lines for that strain
        if (start_pos is not None and start < start_pos) or (end_pos is not None and end > end_pos):
            exclusion_dict[strain] += 1 
            continue

        # Create the event identifier
        event = f"{start}-{end}-{indel_type}"

        # Initialize the strain entry in the dictionary if not already present
        if strain not in presence_absence_dict:
            presence_absence_dict[strain] = {}

        # Mark the presence of the event for the strain
        presence_absence_dict[strain][event] = 1

    # Convert the dictionary to a dataframe
    presence_absence_df = pd.DataFrame.from_dict(presence_absence_dict, orient='index').fillna(0).astype(int)

    # Check if all lines of a given strain were excluded and append the info to the exclusion file
    for key in strain_count:
        if strain_count[key] == exclusion_dict[key]:
            append_line(exclusions, key)


    return presence_absence_df

##################################################################
def append_line(filename, line):
    with open(filename, 'a') as file:
        file.write(line + '\n')

##################################################################
def main():
    parser = argparse.ArgumentParser(description='Concatenate BED files produced by IDIndel.py from a directory, also modifies the first column to match the sample ID.')
    parser.add_argument('-i','--input', metavar='directory', type=str, required=True, help='Input directory containing BED files')
    parser.add_argument('-p','--prefix', metavar='string', type=str, required=True, help='Output prefix for the concatenated file')
    parser.add_argument('-s','--start_pos', metavar='integer', type=int, help='Exclude events that start before this position')
    parser.add_argument('-e','--end_pos', metavar='integer', type=int, help='Exclude events that start after this position')

    args = parser.parse_args()

    no_indel = args.prefix+'_no-Indels.txt'
    with open(no_indel, "w") as ni:
        pass

    indel_out_of_bounds = args.prefix+'_Out-of-bounds-Indels.txt'
    with open(indel_out_of_bounds, "w") as oob:
        pass

    combined_bed = concatenate_bed_files(args.input, no_indel)

    presence_absence_matrix = build_presence_absence_matrix(combined_bed,indel_out_of_bounds, args.start_pos, args.end_pos)
    presence_absence_matrix.to_csv(args.prefix + '_PA_matrix.csv')

    print("Matrix construction completed\n")
    
if __name__ == "__main__":
    main()
