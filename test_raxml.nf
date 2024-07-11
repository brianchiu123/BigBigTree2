#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.file_path="test_folder"
params.placement="$baseDir/scripts/placement.py"



process step4_3_2_placement_trees {
    cpus 4 
    
    input:
    path place_py
    path ph_files_process
    path fasta_files_process
    path unclustered_files
    


    script:
    """
    python3 ${place_py} ${ph_files_process} ${fasta_files_process} ${unclustered_files}
    """
}

workflow{
    place_py = file(params.placement)
    file_list = file("${params.file_path}/*").collect()
    cluster_files = file_list.findAll { it.name.startsWith('cluster') }

    ph_files = cluster_files.findAll { it.name.endsWith('.ph') }
    fasta_files = cluster_files.findAll { it.name.endsWith('.fasta') }

    unclustered_files = file("${params.file_path}/unclustered_seqs.fasta")

    step4_3_2_placement_trees(place_py,ph_files.collect(), fasta_files.collect(), unclustered_files)

}