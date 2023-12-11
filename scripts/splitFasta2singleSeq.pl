#!/usr/bin/perl -w 

use Bio::SeqIO;
my $INPUT=$ARGV[0];
my $COUNT=1;

my $seqs = Bio::SeqIO->new(-file => "<$INPUT");
while (my $seq = $seqs->next_seq) 
{
  my $OUTPUT="tmp$COUNT.fasta";
  my $out = Bio::SeqIO->new(-file => ">$OUTPUT");
  $out->write_seq($seq);
  $COUNT++;
}

exit;
