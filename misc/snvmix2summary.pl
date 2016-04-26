#!/usr/bin/perl
# Script to threshold output values from SNVMix2
# Author: Rodrigo Goya, 05/2009
#
# Copyright (c) 2009, by Sohrab Shah <sshah@bccrc.ca> and Rodrigo Goya <rgoya@bcgsc.ca>
#
# Licensed under the MIT license: http://www.opensource.org/licenses/mit-license.php */

use strict;

use Getopt::Std;
my $opt_string = 'hi:c:t:';
my %opt;
getopts( "$opt_string", \%opt ) or usage();
usage() if $opt{h};
my $SNVMIX_FILE = "-";
$SNVMIX_FILE = $opt{i} if $opt{i};
my $TYPE = 2;
$TYPE = $opt{c} if $opt{c};
my $THRESHOLD = 0;
$THRESHOLD = $opt{t} if $opt{t};
if($TYPE != 2 && $TYPE != 3) { die("ERROR: Unknown class TYPE\n"); }

print STDERR "Reading from ".($SNVMIX_FILE eq "-" ? "STDIN" : $SNVMIX_FILE)."\n";
print STDERR "Calculating for max between AA".($TYPE == 2 ? " and {AB u BB}" : ", AB and BB")."\n";
if($THRESHOLD) {
	print STDERR "Applying threshold of $THRESHOLD, reporting only if ".($TYPE == 2 ? "P{AB u BB}" : "(P{AB} || P{BB})")." >= $THRESHOLD\n";
}

open(INPUT, "<$SNVMIX_FILE") || die("ERROR: Could not open '$SNVMIX_FILE' for reading\n");
while(<INPUT>) {
	chomp;
	s///;
	my $line = $_;
	my ($chr_pos, $ref, $nref, $call_str) = split(/\t/, $line);
	my ($ref_num, $nref_num, $pAA, $pAB, $pBB, $call) = split(/,/, $call_str);
	my $snv = 0;
	if($TYPE == 2) {
		if($pAA < ($pAB + $pBB)) {
			if( ($pAB + $pBB) >= $THRESHOLD) {
				$snv = 1;
			}
		}
	} elsif($TYPE == 3) {
		if($call == 2 || $call == 3) {
			if( $pAB >= $THRESHOLD || $pBB >= $THRESHOLD) {
				$snv = 1;
			}
		}
	} else {
		die("ERROR, and a weird one, script shouldn't even BE in here...\n");
	}
	if($snv) {
		#print "$chr_pos\t$ref\t".( $snv ? $nref : "-")."\t$snv\n";
		print "$line\n";
	}
}
close(INPUT);

sub usage() {
	print "Syntax:\n";
	print "$0 [-i <file>] -c <TYPE> [-t <THRESHOLD>]\n";
	print "\tIf file not given, STDIN is read\n";
	print "\tTYPE is the number of classes to consider\n";
	print "\t\t'2'\tconsiders only AA and {AB U BB} (default)\n";
	print "\t\t'3'\tconsiders AA, AB and BB\n";
	print "\tIf -t THRESHOLD is given, then SNVs will be reported\n";
	print "\twhen the selected probability exceeds this\n";
	exit;
}
