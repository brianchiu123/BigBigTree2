from Bio import SeqIO
import os

# 要合併的FASTA文件列表
fasta_files = ['file1.fasta', 'file2.fasta', 'file3.fasta']  # 請替換為你的文件名

# 使用字典來存儲每個序列ID及其對應的最長序列
sequence_dict = {}

for fasta_file in fasta_files:
    # 確保文件存在
    if not os.path.exists(fasta_file):
        print(f"File {fasta_file} not found.")
        continue
    
    # 讀取FASTA文件
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        seq = str(record.seq)
        
        # 如果這個序列ID已經在字典中，並且現有的序列比新的短，則更新
        if seq_id in sequence_dict:
            if len(sequence_dict[seq_id]) < len(seq):
                sequence_dict[seq_id] = seq
        else:
            sequence_dict[seq_id] = seq

# 創建一個新的序列列表，用於寫入新的FASTA文件
new_sequences = [SeqIO.SeqRecord(SeqIO.Seq(seq), id=seq_id) for seq_id, seq in sequence_dict.items()]

# 寫入新的FASTA文件
with open("merged_sequences.fasta", "w") as output_handle:
    SeqIO.write(new_sequences, output_handle, "fasta")

print("合併完成。")
