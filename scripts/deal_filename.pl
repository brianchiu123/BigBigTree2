#!/usr/bin/perl -w 
use Bio::Seq;
use Bio::SeqIO;
foreach $argnum (0 .. $#ARGV){
$aa = $ARGV[$argnum];



$newa=$aa;
if($aa=~ m/^cluster/)
{$newa=~ s/fasta_aln_aa/aa/ ;}
else {
$newa=~ s/fasta_aln_aa/aln_aa/ ;
}
 rename $aa, $newa ;

}

