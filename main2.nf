nextflow.enable.dsl=2

params.aa="$baseDir/example/aa1.fasta"
params.speciesTree="$baseDir/example/test_species.ph"
params.nn="$baseDir/example/nn1.fasta"
params.cluster_dir='res_dic/cluster'
params.msa_mode='tcoffee'
params.tcoffee_mode='fmcoffee'
params.cluster_Number='5'
params.py_diff="$baseDir/scripts/fasta_dif.py"
params.logfile="$baseDir/nextflow.log"
params.tree_mode="phyml"
params.file_path=".file_path"
params.output="$baseDir/output"


log.info """\
         R N A T O Y   P I P E L I N E    
         =============================
         aa: ${params.aa}
         speciesTree : ${params.speciesTree}
         nn : ${params.nn}
         cluster_dir: ${params.cluster_dir}
		
         """
         .stripIndent()



process step0_check_fasta_diff{
    input:

    file SEQ_aa 
    file SEQ_nn
    file diff_in 

    output:
    stdout emit: result


    /// script:
    """
    python3 $diff_in $SEQ_aa $SEQ_nn
        
    """
}

process step1_1cluster {
	 
    input:	
    file SEQ_fasta_aa 

    output:
    path 'pairwiseSim.fasta' , emit : SIMILARITY_FILE 
    stdout emit : result1_1

    script :

    """
    step1_1  $SEQ_fasta_aa 
	
    """
}

process step1_2dif  {

    publishDir "${params.output}/cluster", pattern: "cluster.txt", mode: 'copy'

    input:
    path SIMILARITY_FILE
    path SEQ_fasta_aa 

    output:
    path 'cluster.txt' , emit :  TEXT_CLUSTER
    path 'clusterTable.csv', emit :  CSV_CLUSTER
    stdout emit: result1_2
	
    """
    step1_2  $SIMILARITY_FILE $SEQ_fasta_aa ${params.cluster_Number}
    """
}

process step2_cluster_to_fasta_aa  {
	
    input:
    path SEQ_fasta_aa 
    path TEXT_CLUSTER 
	 
    output:
    path '*.fasta' ,emit : c
    stdout emit: result2_1	 

    """
    step2_cluster-2-fasta-4-seqFasta.pl $SEQ_fasta_aa $TEXT_CLUSTER
    """
}

process step2_cluster_to_fasta_nn  {
    input:
    path SEQ_fasta_nn 
    path TEXT_CLUSTER 

    output:
    path '*.fasta' ,emit: n
    stdout emit: result2_2	

    """
    step2_cluster-2-fasta-4-seqFasta.pl $SEQ_fasta_nn $TEXT_CLUSTER
    """

}

process step3_1_alignment_aa {
	 		
    input:
    path aa_fasta_channel 
	
    output:
    path '*.fasta_aln_aa' ,emit: step_3_1
    stdout emit: result3_1a

    when:
        params.msa_mode == 'tcoffee'

    script:
    """	
    t_coffee ${aa_fasta_channel} -multi_core no +keep_name -output fasta_aln  -mode ${params.tcoffee_mode}  -outfile  "${aa_fasta_channel}_aln_aa"
	
    """	 

}

process step3_1_alignment_nn {
	
    input:
    path alnaa 
    path nnf 

    output:
    path '*.fasta_aln_nn' ,emit: step3_1_nn
    stdout emit: result3_1n

	when:
    params.msa_mode == 'tcoffee' 

    script:
    """	
    t_coffee -other_pg seq_reformat -in $nnf   -in2 $nnf"_aln_aa"  -action +thread_dna_on_prot_aln -output fasta >   $nnf"_aln_nn"
			
    """	 

}

//mafft
process step3_1_alignment_aa_mafft { 

    label 'mafft'

    input:
    path aa_fasta_channel 
    
    output:
    path '*.fasta_aln_aa' ,emit: step3_1_mafft
    stdout emit: result3_1a_mafft

    when:
        params.msa_mode == 'mafft'

    script:
    """ 
    mafft  --auto --inputorder ${aa_fasta_channel} > ${aa_fasta_channel}_aln_aa    
    """  

}

process step3_1_alignment_nn_mafft {

    label 'mafft'

    input:
    path alnaa 
    path nnf 

    output:
    path '*.fasta_aln_nn' , emit : step3_1_mafft_nn
    stdout emit: result3_1n_mafft

    when:
        params.msa_mode == 'mafft'

/*     
    script:
    """
    mafft --auto --inputorder --anysymbol ${nnf} > ${nnf}_aln_nn 
    """
*/

    script:
    """
    t_coffee -other_pg seq_reformat -in $nnf   -in2 $nnf"_aln_aa"  -action +thread_dna_on_prot_aln -output fasta >   $nnf"_aln_nn"
    """

}

process step3_2_deal_filename{

    input:
    file aa

    // file p 

    output:
    path '*.aln_aa' ,emit: alnfa
    stdout emit:result3_2
    """
    deal_filename.pl $aa
	
    """

}

process step3_2_deal_duplicate{

    input:
    file nn 
    file aa 
    // file p from path_file

    output:
    path "*.fasta_ali_nn" ,emit: aln_nn_tocon
    stdout emit :result3_3
    script:
		
    """
    deal-duplicateID.pl $aa

    """

}

process step3_2_concatenate{

    publishDir "${params.output}/concatenation", pattern: "*_aln", mode: 'copy'

    input:
    path nn_tocon 
    // file p from path_file

    output:
    path 'concatenation.fasta_aln' ,emit: concatenate_output
    stdout emit :result3_4
	
    """
    step3-2_concateAlign  "concatenation.fasta_aln"
    """

}

process step4_1_produce_treebest {

    input:
    path con_fasta_aln 
    //file p from path_file

    output:
    path 'concatenation.ph' ,emit : con_ph_best
    stdout emit : result4_1

    when:
        params.tree_mode == 'treebest'

    script:
    """
    treebest best -o concatenation.ph concatenation.fasta_aln
    """

}

process step4_1_produce_tree_phyml {

    input:
    path con_fasta_aln 
    //file p from path_file

    output:
    path 'concatenation.ph' ,emit: con_ph_phy
    stdout emit:result4_2

    when:
        params.tree_mode == 'phyml'

    script:
    """
    touch tempp.code
    touch tempp
    t_coffee -other_pg seq_reformat -in concatenation.fasta_aln -output code_name>  tempp.code
    t_coffee -other_pg seq_reformat -code  tempp.code -in  concatenation.fasta_aln -output phylip > tempp
    phyml -i tempp -b 0
    postprocess-4-phyml.pl tempp concatenation.ph
    """

}

process step4_1_produce_raxml {
    label 'raxml'
    input:
    path con_fasta_aln 
    //file p from path_file

    output:
    path 'concatenation.ph' ,emit: con_ph_raxml

    when:
        params.tree_mode == 'raxml'

    script:
    """
    raxml-ng --msa concatenation.fasta_aln --model GTR+G --thread 4 --seed 2
    mv *.bestTree concatenation.ph
    """
}

process step4_2_deal_cluster{
	
    input:
    path nnn 
    //file p from path_file

    output:
    path 'cluster*.aln_nn' ,emit :aln_4_2
    stdout emit : result4_3
    """
    deal_clustername.pl $nnn
    """

}


process step4_2_produce_tree {

    errorStrategy 'ignore'

    input:
    path aln_nn 
    path con 
    path speciesTree 
    //file p from path_file

    output:
    path '*.ph' ,emit: clu_ph
    stdout emit: result4_4

    when:
        params.tree_mode == 'treebest'

    script:
    """
    treebest best -o ${aln_nn.getBaseName()}".ph" $aln_nn -f $speciesTree

    """

}

process step4_2_produce_tree_phyml {
    errorStrategy 'ignore'
    input:
    path aln_nn 
    path con 
    path speciesTree
    //file p from path_file

    output:
    path '*.ph' ,emit: clu_ph_phy
    stdout emit: result4_5

    when:
        params.tree_mode == 'phyml'

    script:
    """
    touch tempp.code
    touch tempp
    t_coffee -other_pg seq_reformat -in $aln_nn -output code_name>  tempp.code
    t_coffee -other_pg seq_reformat -code  tempp.code -in  $aln_nn -output phylip > tempp
    phyml -i tempp -b 0
    postprocess-4-phyml.pl tempp ${aln_nn.getBaseName()}.ph

    """

}

process step4_2_produce_tree_raxml {
    //label 'raxml'
    errorStrategy 'ignore'

    input:
    path aln_nn 
    path con 
    path speciesTree 
    //file p from path_file

    output:
    path '*.ph' ,emit:clu_ph_raxml

    when:
        params.tree_mode == 'raxml'

    //use simphy to generate subtree
    script:
    """
    touch tempp.code
    touch tempp
    t_coffee -other_pg seq_reformat -in $aln_nn -output code_name>  tempp.code
    t_coffee -other_pg seq_reformat -code  tempp.code -in  $aln_nn -output phylip > tempp
    phyml -i tempp -b 0
    postprocess-4-phyml.pl tempp ${aln_nn.getBaseName()}.ph
    """

/*
    script:
    """
    treebest best -o ${aln_nn.getBaseName()}".ph" $aln_nn -f $speciesTree
    """
*/

/*
    script:
    """
    raxml-ng --msa $aln_nn --model GTR+G --thread 4 --seed 2
    mv *.bestTree ${aln_nn.getBaseName()}.ph
    """
*/

}

process step4_3_produce_tree {

    publishDir "${params.output}/final_tree", pattern: "*.ph", mode: 'copy'

    input:
    path cluster_ph 
    path con 
    //file p from path_file

    output:
    path 'final.ph' ,emit: final_result
    stdout emit:result4_6 

    script:	
    """
    mergeGroup2bigTree concatenation.ph final.ph 
    """

}


workflow{
    aa_file=file(params.aa)
    speciesTree_file=file(params.speciesTree)
    nn_file=file(params.nn)
    diff= file(params.py_diff)
    log_file=file(params.logfile)
    path_file=file(params.file_path)

    step0_check_fasta_diff(aa_file,nn_file,diff)
    step0_check_fasta_diff.out.result.subscribe {println it}

    step1_1cluster(aa_file)
    SIMILARITY_FILE  = step1_1cluster.out.SIMILARITY_FILE

    step1_2dif(SIMILARITY_FILE,aa_file)
    TEXT_CLUSTER = step1_2dif.out.TEXT_CLUSTER
    CSV_CLUSTER = step1_2dif.out.CSV_CLUSTER



    step2_cluster_to_fasta_aa(aa_file,TEXT_CLUSTER)
    cluster_fasta_aa = step2_cluster_to_fasta_aa.out.c
    cluster_fasta_aa_mafft = step2_cluster_to_fasta_aa.out.c



    step2_cluster_to_fasta_nn(nn_file,TEXT_CLUSTER)
    cluster_fasta_nn = step2_cluster_to_fasta_nn.out.n
    cluster_fasta_nn_mafft = step2_cluster_to_fasta_nn.out.n



    step3_1_alignment_aa(cluster_fasta_aa.flatten())
    aln_fasta_aa = step3_1_alignment_aa.out.step_3_1
    aln_fasta_aa_3_2 = step3_1_alignment_aa.out.step_3_1

    step3_1_alignment_nn(aln_fasta_aa.collect(),cluster_fasta_nn.flatten())
    aln_fasta_nn = step3_1_alignment_nn.out.step3_1_nn
    aln_fasta_nn_4_1 = step3_1_alignment_nn.out.step3_1_nn


    step3_1_alignment_aa_mafft(cluster_fasta_aa.flatten())      
    aln_fasta_aa_mafft = step3_1_alignment_aa_mafft.out.step3_1_mafft
    aln_fasta_aa_3_2_mafft = step3_1_alignment_aa_mafft.out.step3_1_mafft

    step3_1_alignment_nn_mafft(aln_fasta_aa_mafft.collect(),cluster_fasta_nn_mafft.flatten())
    aln_fasta_nn_mafft = step3_1_alignment_nn_mafft.out.step3_1_mafft_nn
    aln_fasta_nn_4_1_mafft = step3_1_alignment_nn_mafft.out.step3_1_mafft_nn

    aln_fasta_aa_3_2_compose = aln_fasta_aa_3_2
        .concat(aln_fasta_aa_3_2_mafft)
    aln_fasta_nn_compose = aln_fasta_nn
        .concat(aln_fasta_nn_mafft)
    aln_fasta_nn_4_1_compose = aln_fasta_nn_4_1
        .concat(aln_fasta_nn_4_1_mafft)

    
    step3_2_deal_filename(aln_fasta_aa_3_2_compose.collect())
    alnfa = step3_2_deal_filename.out.alnfa

    step3_2_deal_duplicate(aln_fasta_nn_compose.collect(),alnfa.flatten())
    aln_nn_tocon = step3_2_deal_duplicate.out.aln_nn_tocon

    step3_2_concatenate(aln_nn_tocon.collect())
    //step3_2_concatenate.out.concatenate_output.view()
    concatenate_aln_best = step3_2_concatenate.out.concatenate_output
    concatenate_aln_phy = step3_2_concatenate.out.concatenate_output
    concatenate_aln_raxml = step3_2_concatenate.out.concatenate_output


    step4_1_produce_treebest(concatenate_aln_best)
    con_ph_best = step4_1_produce_treebest.out.con_ph_best

    // ???
    step4_1_produce_tree_phyml(concatenate_aln_phy)
    con_ph_phy = step4_1_produce_tree_phyml.out.con_ph_phy

    //???
    step4_1_produce_raxml(concatenate_aln_raxml)
    con_ph_raxml = step4_1_produce_raxml.out.con_ph_raxml

    step4_2_deal_cluster(aln_fasta_nn_4_1_compose.collect())
    aln_4_2 = step4_2_deal_cluster.out.aln_4_2
    aln_4_2_phy = step4_2_deal_cluster.out.aln_4_2
    aln_4_2_raxml = step4_2_deal_cluster.out.aln_4_2

    //merge
    con_ph = con_ph_best
    .concat(con_ph_phy,con_ph_raxml)

    step4_2_produce_tree(aln_4_2.flatten(),con_ph_best,speciesTree_file)
    clu_ph = step4_2_produce_tree.out.clu_ph

    // phyml??? 
    step4_2_produce_tree_phyml(aln_4_2_phy.flatten(),con_ph_phy,speciesTree_file)
    clu_ph_phy = step4_2_produce_tree_phyml.out.clu_ph_phy

    // raxml???
    step4_2_produce_tree_raxml(aln_4_2_raxml.flatten(),con_ph_raxml,speciesTree_file)
    clu_ph_raxml = step4_2_produce_tree_raxml.out.clu_ph_raxml


    //merge
    clu_p = clu_ph
    .concat(clu_ph_phy,clu_ph_raxml)

    step4_3_produce_tree(clu_p.collect(),con_ph)

    final_result = step4_3_produce_tree.out.final_result


}

