# 定義一個簡單的映射表，用於將蛋白質的胺基酸轉換為一個對應的DNA三连体（一種可能性）
codon_table = {
    'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT',
    'C': 'TGT', 'Q': 'CAA', 'E': 'GAA', 'G': 'GGT',
    'H': 'CAT', 'I': 'ATT', 'L': 'TTA', 'K': 'AAA',
    'M': 'ATG', 'F': 'TTT', 'P': 'CCT', 'S': 'TCT',
    'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
    '*': 'TAA', # 停止密碼子
    'X': 'NAA' # 用於表示任意或未知胺基酸
}

def protein_to_cds(protein_seq):
    """將蛋白質序列轉換為CDS序列"""
    return ''.join([codon_table.get(aa, 'NNN') for aa in protein_seq])  # 若找不到對應的codon，使用'NNN'

def convert_fasta_protein_to_cds(input_file_path, output_file_path):
    """讀取FASTA格式的蛋白質序列文件，並將其轉換為CDS序列，保存到新文件中"""
    with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
        for line in infile:
            line = line.strip()
            if line.startswith('>'):  # FASTA header
                outfile.write(line + '\n')
            else:
                cds_seq = protein_to_cds(line)
                outfile.write(cds_seq + '\n')

# 定義輸入和輸出檔案路徑
input_file_path = 'tree.PTHR24416_SF263.ENSGT00940000165255.seq.prot.fa' # 這裡填入您的輸入檔案路徑
output_file_path = 'tree.PTHR24416_SF263.ENSGT00940000165255.V2.seq.prot.fa' # 這裡填入您希望保存的輸出檔案路徑

# 執行轉換
convert_fasta_protein_to_cds(input_file_path, output_file_path)

