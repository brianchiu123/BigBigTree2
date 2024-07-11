from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def translate_dna_to_protein(dna_seq):
    return dna_seq.translate()

def find_mismatch(dna_seq, protein_seq):
    translated_protein_seq = translate_dna_to_protein(dna_seq)
    return translated_protein_seq != protein_seq

def fix_mismatch(dna_records, protein_records):
    fixed_dna_records = []
    for protein_record in protein_records:
        protein_seq = str(protein_record.seq)
        original_dna_record = dna_records.get(protein_record.id)
        
        if original_dna_record is None:
            print(f"{protein_record.id}: Not found in original DNA file")
            continue
        
        dna_seq = str(original_dna_record.seq)
        fixed_dna_seq = ""
        dna_index = 0
        
        for aa in protein_seq:
            if aa == "-":
                fixed_dna_seq += "---"
            else:
                codon = dna_seq[dna_index:dna_index+3]
                fixed_dna_seq += codon
                dna_index += 3
        
        fixed_dna_records.append(SeqRecord(Seq(fixed_dna_seq), id=protein_record.id, description=""))
    
    return fixed_dna_records

def restore_original_dna(fixed_dna_records):
    restored_dna_records = []
    for record in fixed_dna_records:
        original_dna_seq = str(record.seq).replace("-", "")
        restored_dna_records.append(SeqRecord(Seq(original_dna_seq), id=record.id, description=""))
    return restored_dna_records


# 读取蛋白质对齐文件
protein_aln_file = "/home/ubuntu/BigBigTree2/work/a2/5dc344447ee732b4818e90d777b049/AccipiterNisus.fasta_aln_aa"
protein_records = list(SeqIO.parse(protein_aln_file, "fasta"))

# 读取原始DNA文件
dna_file = "/home/ubuntu/BigBigTree2/work/a2/5dc344447ee732b4818e90d777b049/AccipiterNisus.fasta"
dna_records = {record.id: record for record in SeqIO.parse(dna_file, "fasta")}

# 修复不一致的DNA序列
fixed_dna_records = fix_mismatch(dna_records, protein_records)

# 保存修复后的DNA文件
fixed_output_file = "AccipiterNisus_fixed.fasta"
SeqIO.write(fixed_dna_records, fixed_output_file, "fasta")
print(f"Fixed DNA sequences have been saved to {fixed_output_file}")

# 恢复原始未对齐的DNA序列
restored_dna_records = restore_original_dna(fixed_dna_records)

# 保存恢复后的DNA文件
restored_output_file = "AccipiterNisus_restored.fasta"
SeqIO.write(restored_dna_records, restored_output_file, "fasta")
print(f"Restored DNA sequences have been saved to {restored_output_file}")