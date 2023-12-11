#!/usr/bin/env perl

$tmp     = $ARGV[0];
$tree    = $ARGV[1];
$code    = $tmp . ".code";
$tmp_tree= $tmp . "_phyml_tree.txt";
$clean   = $tmp . "*";

`cp $tmp_tree $tree`;
open HASH_FILE, "$code" or die $!;
my $tmp=<HASH_FILE>;
foreach $line (<HASH_FILE>)
{
  chomp($line);
  my @tmps = split(" ", $line);
  my $result = `sed -i s/$tmps[1]:/$tmps[0]:/ $tree`;
}
close(HASH_FILE);

`rm $clean`;
