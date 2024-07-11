from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def translate_dna_to_protein(dna_seq):
    return Seq(dna_seq).translate()

def find_codon_for_aa(dna_seq, aa, start_pos):
    while start_pos < len(dna_seq) - 2:
        codon = dna_seq[start_pos:start_pos+3]
        if translate_dna_to_protein(codon) == aa:
            return codon, start_pos
        start_pos += 3
    return None, start_pos

# 读取蛋白质对齐文件
protein_aln_file = "/home/ubuntu/BigBigTree2/work/a2/5dc344447ee732b4818e90d777b049/AccipiterNisus.fasta_aln_aa"
protein_records = list(SeqIO.parse(protein_aln_file, "fasta"))

# 读取原始DNA文件
dna_file = "/home/ubuntu/BigBigTree2/work/a2/5dc344447ee732b4818e90d777b049/AccipiterNisus.fasta"
dna_records = {record.id: record for record in SeqIO.parse(dna_file, "fasta")}

# 反向翻译蛋白质对齐文件
translated_dna_records = []
for protein_record in protein_records:
    protein_seq = protein_record.seq
    dna_seq = ""
    dna_index = 0
    original_dna_record = dna_records.get(protein_record.id)
    if original_dna_record is None:
        print(f"{protein_record.id}: Not found in original DNA file")
        continue
    
    original_dna_seq = str(original_dna_record.seq)
    
    for aa in protein_seq:
        if aa == "-":
            dna_seq += "---"
        else:
            codon, dna_index = find_codon_for_aa(original_dna_seq, aa, dna_index)
            if codon:
                dna_seq += codon
                dna_index += 3
            else:
                print(f"{protein_record.id}: Codon not found for {aa} starting from position {dna_index}")
                break

    translated_dna_records.append(SeqRecord(Seq(dna_seq), id=protein_record.id, description=""))

# 比较反向翻译的DNA与原始DNA
for translated_record in translated_dna_records:
    original_record = dna_records.get(translated_record.id)
    if original_record:
        if str(translated_record.seq).replace("-", "") == str(original_record.seq):
            print(f"{translated_record.id}: Match")
        else:
            print(f"{translated_record.id}: Mismatch")
    else:
        print(f"{translated_record.id}: Not found in original DNA file")
