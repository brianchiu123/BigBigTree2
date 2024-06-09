
import csv

def read_tsv(filename):
    with open(filename, 'r', encoding='utf-8') as file:
        reader = csv.reader(file, delimiter='\t')
        rows = [row for row in reader]

    return rows



def read_fasta(filename,transform_data):
    with open(filename, 'r') as file:
        sequences = {}
        sequence_name = None
        sequence_data = ""
        transformIdDate = []
        tempID = ""
        
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
                        tempID = sequence_name + "_" + row[2].replace("_","")
            else:
                sequence_data += line
                
        
        # Handle the last sequence
        if sequence_name:
            transformIdDate.append([tempID,sequence_data])
            sequences[sequence_name] = sequence_data

    return transformIdDate

def write_fasta(data, filename):
    with open(filename, 'w') as file:
        for item in data:
            sequence_name, sequence = item
            file.write(f">{sequence_name}\n{sequence}\n")




tsv_filename = 'INSR_testData/tree.PTHR24416_SF263.ENSGT00940000165255.gene_list.tsv'

transform_data = read_tsv(tsv_filename )

aa_file = "INSR_testData/tree.PTHR24416_SF263.ENSGT00940000165255.seq.prot.fa"
aa_sequences = read_fasta(aa_file,transform_data)
aa_output = "INSR_transformData/tree.PTHR24416_SF263.ENSGT00940000165255.seq.prot.fa"
write_fasta(aa_sequences,aa_output)

nn_file = "INSR_testData/tree.PTHR24416_SF263.ENSGT00940000165255.seq.cds.fa"
nn_sequences = read_fasta(nn_file,transform_data)
nn_output = "INSR_transformData/tree.PTHR24416_SF263.ENSGT00940000165255.seq.cds.fa"
write_fasta(nn_sequences,nn_output)





