#!/usr/bin/perl -w 
use Bio::Seq;
use Bio::SeqIO;
foreach $argnum (0 .. $#ARGV){
$aa = $ARGV[$argnum];



$newa=$aa;
if($aa=~ m/^cluster/)
{$newa=~ s/fasta_aln_nn/aln_nn/ ;}
else {
$newa=~ s/fasta_aln_nn/nn/ ;
}
 rename $aa, $newa ;

}

