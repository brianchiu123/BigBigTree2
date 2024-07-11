import sys
import subprocess
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

def fix_mismatch_and_restore(dna_file, protein_aln_file, fixed_dna_file):
    protein_records = list(SeqIO.parse(protein_aln_file, "fasta"))
    dna_records = {record.id: record for record in SeqIO.parse(dna_file, "fasta")}
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
                    print(f"{protein_record.id}: Codon not found for {aa} starting from position {dna_index}. Using 'NNN' as placeholder.")
                    dna_seq += "NNN"
                    dna_index += 3

        translated_dna_records.append(SeqRecord(Seq(dna_seq), id=protein_record.id, description=protein_record.description))

    # Restore original DNA sequences by removing gaps
    restored_dna_records = []
    for record in translated_dna_records:
        original_dna_seq = str(record.seq).replace("-", "")
        restored_dna_records.append(SeqRecord(Seq(original_dna_seq), id=record.id, description=record.description))

    # 保存修正后的 DNA 文件
    SeqIO.write(translated_dna_records, fixed_dna_file, "fasta")
    # 保存恢复的原始 DNA 文件
    restored_dna_file = fixed_dna_file.replace("fixed_", "restored_")
    SeqIO.write(restored_dna_records, restored_dna_file, "fasta")

    return restored_dna_file

def run_t_coffee(fixed_dna_file, protein_aln_file, output_file):
    cmd = f"t_coffee -other_pg seq_reformat -in {fixed_dna_file} -in2 {protein_aln_file} -action +thread_dna_on_prot_aln -output fasta > {output_file}"
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        raise RuntimeError(f"T-Coffee failed for {fixed_dna_file} with {protein_aln_file}")

if __name__ == "__main__":
    dna_file = sys.argv[1]
    protein_aln_file = sys.argv[2]
    fixed_dna_file = sys.argv[3]
    output_file = sys.argv[4]

    restored_dna_file = fix_mismatch_and_restore(dna_file, protein_aln_file, fixed_dna_file)
    run_t_coffee(restored_dna_file, protein_aln_file, output_file)
