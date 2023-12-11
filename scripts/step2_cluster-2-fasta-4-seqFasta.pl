#!/usr/bin/perl -w 
use Bio::Seq;
use Bio::SeqIO;

my %species_hash= ();
my %cluster_hash = ();
my %count_hash = ();

if (@ARGV==0){
        print "Usage: cluster-2-fasta.pl sequence_file cluster_file\n";
        exit (1);  
}
my $seqs	= $ARGV[0];
my $cluster	= $ARGV[1];


print "\tseq     file= $seqs\n";
print "\tcluster file= $cluster\n";

my $inseq = Bio::SeqIO->new(-file => "<$seqs");
while (my $seq = $inseq->next_seq) {
	my @tmps=split('_', $seq->id);
	$species=$tmps[$#tmps];
	
#get cluster ID from Cluster file
#for a species, more than one sequences,ie. seqX and seqY,belong to one cluster,ie. cluster1, naming rule:
#cluster1_1
#cluster1_2
	$query_id=$seq->id;
	$cluster_id = `grep -n $query_id $cluster|awk -F: '{printf \"\%d\", \$1}'`;
	$check_index="$species.$cluster_id";
	$count_hash{$check_index}++;
	$outcluster_id="cluster".$cluster_id."_".$count_hash{$check_index};
			
        $obj4species = Bio::Seq->new(-display_id => $outcluster_id, -desc => $seq->id, -seq => $seq->seq);
        $obj4cluster = Bio::Seq->new(-display_id => $seq->id, -seq => $seq->seq);
		
	push @{$species_hash{$species}}, $obj4species;
	push @{$cluster_hash{$cluster_id}}, $obj4cluster;	
}

#output fasta file by species, sort by group ID
for $species ( keys %species_hash ) {
	@{$species_hash{$species}} = sort { $a->id cmp $b->id } @{$species_hash{$species}}; 
}
for $species ( keys %species_hash ) {
	$out = Bio::SeqIO->new(-file => ">$species.fasta" , '-format' => 'fasta');
	foreach (@{$species_hash{$species}}) {
		$out->write_seq($_);
    	}
}

#output fasta file by group
for $group ( keys %cluster_hash ) {
	$out = Bio::SeqIO->new(-file => ">cluster$group.fasta" , '-format' => 'fasta');
	foreach (@{$cluster_hash{$group}}) {
		$out->write_seq($_);
    	}
}

exit(0);
