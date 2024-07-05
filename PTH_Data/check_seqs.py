from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Function to translate a DNA sequence to a protein sequence
def translate_dna_to_protein(dna_seq):
    return dna_seq.translate(to_stop=True)

# Function to read sequences from a FASTA file
def read_fasta(file_path):
    records = SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))
    return records

# Function to check if protein sequences match and correct mismatches
def check_and_correct_sequences(gene_sequences, protein_sequences):
    corrected_sequences = []
    for gene_id, gene_seq in gene_sequences.items():
        protein_seq = translate_dna_to_protein(gene_seq.seq)
        if gene_id in protein_sequences:
            if protein_seq == protein_sequences[gene_id].seq:
                corrected_sequences.append(gene_seq)
            else:
                print(f"Mismatch found for {gene_id}")
                print(f"Translated protein: {protein_seq}")
                print(f"Given protein: {protein_sequences[gene_id].seq}")
                # Correcting the sequence by using the translated protein sequence
                corrected_seq = SeqRecord(protein_sequences[gene_id].seq.back_transcribe(), id=gene_id, description="corrected")
                corrected_sequences.append(corrected_seq)
        else:
            print(f"No corresponding protein sequence for {gene_id}")
            corrected_sequences.append(gene_seq)
    return corrected_sequences

# Usage example
gene_fasta_file = "PTH1000/PTH1000_nn_sequences.fa"
protein_fasta_file = "PTH1000/PTH1000_aa_sequences.fa"

# Read sequences from FASTA files
gene_sequences = read_fasta(gene_fasta_file)
protein_sequences = read_fasta(protein_fasta_file)

# Check and correct sequences
corrected_sequences = check_and_correct_sequences(gene_sequences, protein_sequences)

# Save corrected sequences to a new FASTA file
with open("corrected_sequences.fasta", "w") as output_handle:
    SeqIO.write(corrected_sequences, output_handle, "fasta")

