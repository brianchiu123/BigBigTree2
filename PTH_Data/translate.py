from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import BiopythonWarning
import warnings

# 忽略Biopython警告
warnings.simplefilter('ignore', BiopythonWarning)

def translate_dna_to_protein(dna_seq):
    # Adjust the length of the sequence to be a multiple of three by adding 'N' if necessary
    if len(dna_seq) % 3 != 0:
        # Add 'N' to make the sequence length a multiple of three
        dna_seq += 'N' * (3 - len(dna_seq) % 3)
    # Translate the DNA sequence to a protein sequence
    return dna_seq.translate(to_stop=True)

def check_and_remove_sequences(dna_fasta, output_file, corrected_dna_fasta, corrected_protein_fasta):
    dna_records = list(SeqIO.parse(dna_fasta, "fasta"))

    corrected_dna_records = []
    corrected_protein_records = []

    with open(output_file, "w") as out_file:
        for dna_record in dna_records:
            dna_id = dna_record.id
            dna_seq = dna_record.seq

            if len(dna_seq) % 3 != 0:
                out_file.write(f"Sequence {dna_id} length is not a multiple of three, added 'N':\n")
                out_file.write(f"Original length: {len(dna_seq)}\n")
                out_file.write(f"Modified length: {len(dna_seq) + (3 - len(dna_seq) % 3)}\n\n")

            translated_seq = translate_dna_to_protein(dna_seq)

            if 'X' in translated_seq:
                out_file.write(f"Sequence {dna_id} contains stop codon:\n")
                out_file.write(f"Translated sequence: {translated_seq}\n\n")
            else:
                corrected_dna_records.append(dna_record)
                corrected_protein_records.append(SeqRecord(translated_seq, id=dna_id, description=""))

    # 將修正後的DNA序列寫入新的FASTA檔案
    SeqIO.write(corrected_dna_records, corrected_dna_fasta, "fasta")
    # 將修正後的蛋白質序列寫入新的FASTA檔案
    SeqIO.write(corrected_protein_records, corrected_protein_fasta, "fasta")

# 使用範例
dna_fasta = "PTH1000/PTH1000_nn_sequences.fa"
output_file = "PTH1000/sequence_issues.txt"

corrected_dna_fasta = "PTH1000/PTH1000_tnn_sequences.fa"
corrected_protein_fasta = "PTH1000/PTH1000_taa_sequences.fa"

check_and_remove_sequences(dna_fasta, output_file, corrected_dna_fasta, corrected_protein_fasta)

