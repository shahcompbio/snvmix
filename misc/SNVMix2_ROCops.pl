#!/usr/bin/perl
# ROC data operations
# Author: Rodrigo Goya, 11/2009
# Copyright (c) 2009, by Sohrab Shah <sshah@bccrc.ca> and Rodrigo Goya <rgoya@bcgsc.ca>
# Licensed under the MIT license: http://www.opensource.org/licenses/mit-license.php */
#
# Based on Tom Fawcetts paper:
# 	An introduction to ROC analysis 
# 	Pattern Recognition Letters 27 (2006) 861–874

use strict;
use Getopt::Std;
my $opt_string = 'hi:o:f:';
my %opt;
getopts( "$opt_string", \%opt ) or usage();
usage() if $opt{h};
my $input = ($opt{i} ? $opt{i} : '-');
my $output = ($opt{o} ? $opt{o} : '-');
my $FDR = ($opt{f} ? $opt{f} : 0.01);


open(INPUT, "<$input") || die("ERROR: could not open input file ($input)\n");
open(OUTPUT, ">$output") || die("ERROR: could not open output file ($output)\n");

my %DATA;
my $i = 0;
my $P = 0;
my $N = 0;
while(<INPUT>) {
	chomp;
	s///;
	my ($val, $class) = split(/\t/, $_);
	if($class == 1) { $P++; }
	elsif($class == 0) { $N++; }
	else { die("ERROR: class should be either 0 or 1, not '$class'\n"); }
	$DATA{$i}{val} = $val;
	$DATA{$i}{class} = $class;
	$i++;
}
close(INPUT);

print STDERR "P = $P\nN = $N\n";

my @L;
my @f;
foreach my $id (sort { $DATA{$b}{val} <=> $DATA{$a}{val} } keys %DATA) {
	push @L, $DATA{$id}{class};
	push @f, $DATA{$id}{val};
	
}
undef(%DATA);

my @R = get_ROC(\@L, \@f);
print STDERR "Roc columns: thr\tfdr\ttpr\n";
for(my $i = 0; $i < scalar(@R); $i++) {
	print OUTPUT "$R[$i]{thr}\t$R[$i]{fdr}\t$R[$i]{tpr}\n";
}

close(OUTPUT);

my $AUC = get_AUC(\@L, \@f);
print STDERR "AUC = $AUC\n";

my $THR = get_THR(\@R, $FDR);
print STDERR "At FDR = $FDR, THR = $THR\n";


#############
# FUNCTIONS #
#############

sub get_ROC {
	my @L = @{shift @_};
	my @f = @{shift @_};

	my $FP = 0;
	my $TP = 0;
	my @R;
	my $fprev = -100000000000000000000;
	for(my $i = 0; $i < scalar(@L); $i++) {
		if($f[$i] != $fprev) {
			push @R, {thr => $f[$i], fdr => $FP/$N, tpr => $TP/$P};
			$fprev = $f[$i];
		}
		if($L[$i]) {
			$TP++;
		} else {
			$FP++;
		}
	}
	push @R, {thr => $f[-1], fdr => $FP/$N, tpr => $TP/$P};
	return @R;
}

sub get_AUC {
	my @L = @{shift @_};
	my @f = @{shift @_};
	my $N = 0;
	my $P = 0;
	my $FP = 0;
	my $TP = 0;
	my $FP_prev = 0;
	my $TP_prev = 0;
	my $AUC = 0;
	my $fprev = -1000000000000000;
	for(my $i = 0; $i < scalar(@L); $i++) {
		if($f[$i] != $fprev) {
			$AUC += trapezoid_area($FP, $FP_prev, $TP, $TP_prev);
			$fprev = $f[$i];
			$FP_prev = $FP;
			$TP_prev = $TP;
		}
		if($L[$i]) {
			$TP++;
			$P++;
		} else {
			$FP++;
			$N++;
		}
	}
	$AUC += trapezoid_area($FP, $FP_prev, $TP, $TP_prev);
	$AUC = $AUC/($N * $P);
	return $AUC;
}

sub trapezoid_area {
	my $base = abs($_[0] - $_[1]);
	my $height_avg = ($_[2] + $_[3]) /2;
	return ($base*$height_avg);
}


sub get_THR {
	my @R = @{shift @_};
	my $FDR = shift @_;
	$i = 0;
	while($R[$i]{fdr} <= $FDR) { $i++; }
	return $R[$i-1]{thr};
}

sub usage() {
	print STDERR "Syntax: $0 -f <FDR> -i <input file> -o <output file>\n";
	exit();
}


