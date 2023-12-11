#!/usr/bin/perl -w 
my $TCOFFEE_CMD = "t_coffee";

$file=$ARGV[0];

  my %hash_2_name=();
  my %hash_2_value=();
  $nn_aln = $file;
  $nn_aln =~ s/aln_aa/fasta_aln_nn/;
 
 
	  print "PROCESS aa_fasta_aln=$file\n";
	  print "        nn_fasta_aln=$nn_aln\n";
	  print "OUTPUT  fasta_aln=$nn_aln\n\n";
		  
	#calculate pairwise sequence similarity   
	@lines=`$TCOFFEE_CMD -other_pg seq_reformat -in $file -output sim|grep AVG`;
	
	#pick up the most similar pairwise sequence similarity
	   foreach $line (@lines){
		   @tmp=split(' ', $line);
		   if($tmp[2] =~ m/(cluster\d+)/)
		   {
			   if ((! $hash_2_name{$1}) || ($tmp[4] > $hash_2_value{$1}))
			   {
				   $hash_2_name{$1}=$tmp[2];
				   $hash_2_value{$1}=$tmp[4];
			   }
		   }
	   }
	
	#extract sequence by similarity   
	   my $list=();
	   foreach $key (keys %hash_2_name)
	   {
		   $list .= "'$hash_2_name{$key}' "; 
	   }
	   system("$TCOFFEE_CMD -other_pg seq_reformat -in $nn_aln -action +keep_name +extract_seq_list $list +rm_gap > tmp.fasta_aln");
	   $nn_aln=~ s/fasta_aln_nn/fasta_ali_nn/;
	   system("sed -i 's/\\(>cluster[0-9]*\\)_[0-9]*/\\1/g' tmp.fasta_aln");
	   system("mv tmp.fasta_aln $nn_aln");  
	 
  

