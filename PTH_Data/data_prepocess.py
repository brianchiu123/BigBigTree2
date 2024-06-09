
import sys
import csv
import pandas as pd
import os


def create_transformed_fasta_files(input_numbers, data_file_path):


    # Load the data file
    data = pd.read_csv(data_file_path, sep='\t')

    def find_name_stable_id_combinations(numbers, data):
        output = []
        for number in numbers.split(","):
            search_term = f"SF{number}"
            matched_rows = data[data['#name'].str.contains(search_term)]
            for _, row in matched_rows.iterrows():
                combined_string = f"{row['#name']}.{row['stable_id']}"
                output.append(combined_string)
        return output

    def read_tsv(filename):
        with open(filename, 'r', encoding='utf-8') as file:
            reader = csv.reader(file, delimiter='\t')
            rows = [row for row in reader]
        return rows
    
    def capitalize_words(s):
        return ''.join(word.capitalize() for word in s.split('_'))

    def read_fasta(filename, transform_data):
        sequences = {}
        sequence_name = None
        sequence_data = ""
        transformIdDate = []
        tempID = ""

        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    if sequence_name:
                        transformIdDate.append([tempID,sequence_data])
                        sequences[sequence_name] = sequence_data
                        sequence_data = ""
                    sequence_name = line[1:]
                    for row in transform_data:
                        if sequence_name == row[0]:
                            capitalized_row = capitalize_words(row[2])
                            tempID = sequence_name + "_" + capitalized_row
                else:
                    sequence_data += line
            
            if sequence_name:  # Handle the last sequence
                transformIdDate.append([tempID,sequence_data])
                sequences[sequence_name] = sequence_data

        return transformIdDate

    # Generate file names based on the name and stable_id combinations
    name_stable_id_combinations = find_name_stable_id_combinations(input_numbers, data)
    combined_aa_sequences = []
    combined_nn_sequences = []

    for item in name_stable_id_combinations:
        print(f'Combine : {item}')
        tsv_filename = f'./INSR-ENSG00000171105-e101/tree.{item}.gene_list.tsv'
        transform_data = read_tsv(tsv_filename)

        # Process amino acid (aa) file
        aa_file = f'./INSR-ENSG00000171105-e101/tree.{item}.seq.prot.fa'
        aa_sequences = read_fasta(aa_file, transform_data)
        combined_aa_sequences.extend(aa_sequences)

        # Process nucleotide (nn) file
        nn_file = f'./INSR-ENSG00000171105-e101/tree.{item}.seq.cds.fa'
        nn_sequences = read_fasta(nn_file, transform_data)
        combined_nn_sequences.extend(nn_sequences)

    # Function to write combined fasta files
    def write_fasta(data, directory, filename):
        os.makedirs(directory, exist_ok=True)
        full_path = os.path.join(directory, filename)
        with open(full_path, 'w') as file:
            for item in data:
                sequence_name, sequence = item
                file.write(f">{sequence_name}\n{sequence}\n")
    output_directory = '/Users/chiuhsienan/Documents/BigBigTree_data/combined_file'
    combined_aa_file = 'combined_aa_sequences.fa'
    combined_nn_file = 'combined_nn_sequences.fa'
    
    write_fasta(combined_aa_sequences, output_directory,combined_aa_file)
    write_fasta(combined_nn_sequences, output_directory,combined_nn_file)

    return os.path.join(output_directory, combined_aa_file), os.path.join(output_directory, combined_nn_file)


# Example usage: 
# combined_aa_file, combined_nn_file = create_transformed_fasta_files("249 162", "./INSR-ENSG00000171105-e101/supertree.PTHR24416.tree_counts.tsv")
# This function will return paths to the combined fasta files.

def count_species_and_ids(fasta_file):
    species_names = []
    ids = []

    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith(">"):
                # Assuming the format is ">id species name"
                parts = line[1:].strip().split("_")
                ids.append(parts[0])
                if parts[1:] not in species_names:
                    species_names.append(parts[1:])

    return len(species_names), len(ids)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_numbers = sys.argv[1]
        supertree_filename = "./INSR-ENSG00000171105-e101/supertree.PTHR24416.tree_counts.tsv"
        combined_aa_file, combined_nn_file = create_transformed_fasta_files(input_numbers, supertree_filename)

        # Count species and IDs in the combined_aa_sequences.fa file
        species_count, id_count = count_species_and_ids(combined_aa_file)
        print(f"Total number of unique species: {species_count}")
        print(f"Total number of unique IDs: {id_count}")
    else:
        print("Usage: python3 data_prepocess.py <numbers>")