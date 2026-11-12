#!/usr/bin/perl
#Used for get SOMATIC varians from VarScan paired VCF
#Retrieve possible somatic variants from normal sample due to seq error. 


use warnings;
use strict;

my $file = $ARGV[0];

open (IN, "<$file");
while(<IN>){
	chomp;
	if ( $_ =~/#/ ){
		print $_."\n";
	} else {
		@_ = split/\t/, $_;
		my $tumor = $_[10]; ####varscan: tumor col.
		my $normal = $_[9]; ####varscan: normal col.
		my @tumor = split/:/, $tumor;
		$tumor[5] =~s/%//;
        $tumor[5] =$tumor[5]/100;
		my @DP4 = split/,/ ,$tumor[-1];
		my @normal = split/:/, $normal;
		$normal[5] =~s/%//;
		$normal[5] =$normal[5]/100;
		## same function with corr_Varsvan_SS-TAG.py ?
		#if ($normal[5] <= 0.01 && $_[7] =~/SS=1/) {
		#	print $_."\n";
		#	next;
		#}
		next if $tumor[4]< 3; ####AD>=3
		next if ( $DP4[2] < 1 || $DP4[3] < 1 ); ####read2 from each strand
		next if ( $DP4[2]/$tumor[4] >= 0.95 || $DP4[2]/$tumor[4] <= 0.05); ####strand bias;ZL: ratio=0.95
		if ( $_[7] =~/SS=2/ ){
			print $_."\n";
		} elsif ( $_[7] =~/SS=1/ ){
			print $_."\n" if ($normal[4] < 2) ||  ( $normal[4] >=2 && $normal[4] < 4 && $tumor[5]/$normal[5] >=5 ); 
		}
	}
}
close IN;

