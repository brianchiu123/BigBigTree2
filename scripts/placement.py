import os
import sys
import json
import subprocess
import re

def run_mafft(unclustered_seq, aligned_cluster, output_file):
    cmd = f"sudo mafft --add {unclustered_seq} --reorder {aligned_cluster} > {output_file}"
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        raise RuntimeError(f"MAFFT failed for {unclustered_seq} with {aligned_cluster}")

def run_raxml(combined_fasta, cluster_tree, output_prefix):
    cmd = f"raxmlHPC-PTHREADS-AVX -T 1 -f v -s {combined_fasta} -t {cluster_tree} -m PROTGAMMAAUTO -n {output_prefix}"
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        raise RuntimeError(f"RAxML failed for {combined_fasta} with {cluster_tree}")
    labelled_tree = f"RAxML_labelledTree.{output_prefix}"
    portable_tree = f"RAxML_portableTree.{output_prefix}.jplace"
    return labelled_tree, portable_tree

def parse_raxml_score(labelled_tree):
    with open(labelled_tree, 'r') as file:
        for line in file:
            if line.startswith('('):
                score_line = line
                break
    score = float(re.sub(r'\[.*?\]', '', score_line).split(':')[-1].strip("();"))
    return score

def parse_jplace_score(jplace_file):
    with open(jplace_file, 'r') as file:
        data = json.load(file)
        best_score = float('-inf')
        for placement in data['placements']:
            for p in placement['p']:
                score = p[1]  # Extracting the likelihood score
                if score > best_score:
                    best_score = score
    return best_score

def filter_clusters_by_size(fasta_files, threshold=3):
    valid_clusters = []
    for fasta_file in fasta_files:
        with open(fasta_file, 'r') as f:
            seq_count = sum(1 for line in f if line.startswith('>'))
            if seq_count > threshold:
                valid_clusters.append(fasta_file)
            else:
                cluster_id = os.path.basename(fasta_file).split('.')[0]
                print(f"Cluster {cluster_id} has {seq_count} sequences, which is <= {threshold}, skipping it.")
    return valid_clusters

def process_labelled_tree(labelled_tree):
    with open(labelled_tree, 'r') as file:
        content = file.read()
    # Remove [ ] and the contents within them
    content = re.sub(r'\[.*?\]', '', content)
    # Remove QUERY___ prefix
    content = re.sub(r'QUERY___', '', content)
    with open(labelled_tree, 'w') as file:
        file.write(content)

def main(ph_files_process, fasta_files_process, unclustered_files):
    # Filter out clusters with less than or equal to 3 sequences
    fasta_files_process = filter_clusters_by_size(fasta_files_process, threshold=3)
    print(f"Filtered FASTA files: {fasta_files_process}")
    
    unclustered_seqs = []
    with open(unclustered_files, 'r') as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    unclustered_seqs.append(seq)
                seq = line
            else:
                seq += line
        if seq:
            unclustered_seqs.append(seq)

    round_num = 1
    processed_unclustered_seqs = []

    while unclustered_seqs:
        mafft_results_dir = f"mafft_results_round_{round_num}"
        os.makedirs(mafft_results_dir, exist_ok=True)
        scores = []
        for unclustered_seq in unclustered_seqs:
            best_score = float('-inf')
            best_tree = ""
            best_fasta = ""
            for fasta_file in fasta_files_process:
                cluster_id = os.path.basename(fasta_file).split('.')[0]
                output_file = f"{mafft_results_dir}/uncluster_temp_{cluster_id}.fasta"
                
                # Write the unclustered_seq to a temp file for MAFFT
                unclustered_seq_file = f"temp_unclustered_seq_{round_num}.fasta"
                with open(unclustered_seq_file, 'w') as temp_f:
                    temp_f.write(unclustered_seq)
                
                run_mafft(unclustered_seq_file, fasta_file, output_file)
                
                combined_fasta = output_file
                cluster_tree = f"{cluster_id}.ph"
                raxml_output_prefix = f"raxml_uncluster_temp_{round_num}_{cluster_id}"
                
                try:
                    labelled_tree, portable_tree = run_raxml(combined_fasta, cluster_tree, raxml_output_prefix)
                    process_labelled_tree(labelled_tree)  # Process the labelled tree to remove []
                    score = parse_jplace_score(portable_tree)  # Updated to parse from jplace file
                    if score > best_score:
                        best_score = score
                        best_tree = labelled_tree
                        best_fasta = combined_fasta
                except RuntimeError as e:
                    print(e)
                    continue
            
            scores.append((best_score, best_tree, best_fasta, unclustered_seq))
        
        # Sort scores and choose the best for each cluster
        scores.sort(reverse=True, key=lambda x: x[0])
        chosen_clusters = set()
        for score, tree, fasta, unclustered_seq in scores:
            cluster_id = os.path.basename(tree).split('_')[-1]
            if cluster_id not in chosen_clusters:
                chosen_clusters.add(cluster_id)
                print(f"Round {round_num}: Best score for cluster {cluster_id} is {score}, with tree {tree} and fasta {fasta}")
                # Replace the original cluster files with the new ones
                os.replace(fasta, f"{cluster_id}.fasta")
                os.replace(tree, f"{cluster_id}.ph")
                processed_unclustered_seqs.append(unclustered_seq)
        
        # Remove processed unclustered sequences
        unclustered_seqs = [seq for seq in unclustered_seqs if seq not in processed_unclustered_seqs]
        
        # Clean up MAFFT results directory
        for file in os.listdir(mafft_results_dir):
            os.remove(os.path.join(mafft_results_dir, file))
        
        round_num += 1
        if not unclustered_seqs:
            break

if __name__ == "__main__":
    ph_files_process = [file for file in sys.argv[1:] if file.endswith('.ph')]
    fasta_files_process = [file for file in sys.argv[1:] if file.endswith('.fasta') and 'uncluster' not in file]
    unclustered_files = sys.argv[-1]
    
    print(f"sys.argv: {sys.argv}")
    print(f"Received files:")
    print(f"PH files: {ph_files_process}")
    print(f"FASTA files: {fasta_files_process}")
    print(f"Unclustered file: {unclustered_files}")
    
    main(ph_files_process, fasta_files_process, unclustered_files)
