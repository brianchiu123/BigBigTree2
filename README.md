# bigbigtree

bigbigtree is a  divide and concatenate strategy for the phylogenetic reconstruction of large orthologous datasets




Compile/Installation 
--------------------

Install Nextflow

	curl -s https://get.nextflow.io | bash

Clone the git repository on your computer with the following command: 

	git clone https://github.com/jmchanglab/bigbigtree.git
    
or

	./nextflow clone jmchanglab/bigbigtree

If you get some permission error:

	chmod -R 777 bigbigtree/


Make sure you have installed python and perl:bio, or you can use docker container below.




Usage 
--------------------

	
	
	compulsory:  
        	--aa	 input of Proteins sequences in FASTA format,the each sequence ID must end with a suffix specifying the Species information separated by a underscore char.For Example:    >XXX_SPECIES there 'XXX' is the sequence ID and 'SPECIES' is the species code.
		
        	--nn 	 input of Nucleotides sequences in FASTA format,the each sequence ID must end with a suffix specifying the Species information separated by a underscore char.For Example:    >XXX_SPECIES there 'XXX' is the sequence ID and 'SPECIES' is the species code.
		
		
		
	optional:	
		--speciesTree		input of species Tree
		--msa_mode		can choose tcoffee,mafft(default:tcoffee)
		--tcoffee_mode		can choose fmcoffee,mcoffee,tcoffee(default:fmcoffee)
		--tree_mode		can choose treebest,phyml(default:treebest)
		--cluster_Number	recommended set to the number of species or larger(default:9)
		--logfile 	
	
        	
For instance:

	nextflow run main.nf --aa 'example/Or_aa.fasta' --nn 'example/Or_nn_v2.fasta' --speciesTree 'example/speciesTree.ph' --msa_mode 'mafft' --tree_mode 'phyml' --logfile '$baseDir/nextflow.log'

	nextflow run main.nf --aa 'experiment/preprocessing/treedepth_0.5/5x8_treedepth_0.5/5x8_aa.fasta' --nn 'experiment/preprocessing/treedepth_0.5/5x8_treedepth_0.5/5x8_nn.fasta' --speciesTree 'experiment/preprocessing/treedepth_0.5/5x8_treedepth_0.5/s_tree.trees' --mas_mode 'mafft' --tree_mode 'raxml'

Docker 
--------------------
If you have installed Docker, you can execute nextflow code with my docker container.Using the command below to get the docker image: 

	docker pull tsaihanlung/mybiocontainer-python-bioperl

Execute nextflow with docker:

	nextflow run main.nf -with-docker
