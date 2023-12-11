#!/usr/bin/perl -w 

use List::MoreUtils qw/ uniq /;

my %GROUP_TABLE;
my @SPECIES_ARR;

my $GROUP_COUNT = 1;
my $filename = $ARGV[0];

open (INPUT, "$filename") || die "cannot read file: $filename!\n";
while (<INPUT>)
{
	@tmp=split(" ",$_);
	@tmp_arrs=split(",",$tmp[6]);
	foreach $seq (@tmp_arrs)
	{
	  if($seq =~ m/_([^\_]+)$/)
	  {
		push(@SPECIES_ARR, $1);
		$GROUP_TABLE{$GROUP_COUNT}{$1}++;
	  }
	}
	$GROUP_COUNT++;
}
close INPUT;

@SPECIES_ARR = uniq @SPECIES_ARR;
@SPECIES_ARR = sort(@SPECIES_ARR);

print "group\\species";
foreach ( @SPECIES_ARR ) {
    print ",$_";
}
print "\n";

for ($i = 1; $i < $GROUP_COUNT; $i++){
  print "$i";
  foreach ( @SPECIES_ARR ) {
    if(exists $GROUP_TABLE{$i}{$_}){
      print ",$GROUP_TABLE{$i}{$_}"
    }
    else{
      print ","
    }
  }
  print "\n";
}
