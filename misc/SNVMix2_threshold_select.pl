#!/usr/bin/perl
# Scripts automates threshold selection for a combination of parameters for SNVMix2
# Author: Rodrigo Goya, 01/2010
#
# Copyright (c) 2009, by Sohrab Shah <sshah@bccrc.ca> and Rodrigo Goya <rgoya@bcgsc.ca>
#
# Licensed under the MIT license: http://www.opensource.org/licenses/mit-license.php */

# Script will:
#   For each combination of model, baseQ and mapQ:
#       Train on <train set>
#       Calculate threshold on <train set> for <desired FDR>
#       Calculate sens, prec and f-measure on <test set> using selected threshold
#       Output lines such as in Table 2
#   Highlight which model/baseQ/mapQ with what threshold gave best f-measure

use strict;
use IPC::Open2;

my $SNVMIX = "SNVMix2";
my $ROCOPS = "SNVMix2_ROCops.pl";
my %VALID_MODELS = ( "mb" => 3, "m" => 2, "b" => 1, "M" => 2, "Mb" => 2, "MB" => 3, "SNVMix1" => 3);
# 0x1 = vary base quality, 0x2 = vary map quality

use Getopt::Std;
my $opt_string = 'hp:t:T:b:m:M:f:w:o:';
my %opt;
getopts( "$opt_string", \%opt ) or usage();
usage() if $opt{h};
my $pileup = ($opt{p} ? $opt{p} : usage());
my $train_set_file = ($opt{t} ? $opt{t} : usage());
my $test_set_file = ($opt{T} ? $opt{T} : usage());
my $opt_bQs = ($opt{b} ? $opt{b} : 0);
my $opt_mQs = ($opt{m} ? $opt{m} : 0);
my $opt_model = ($opt{M} ? $opt{M} : usage());
my $FDR = ($opt{f} ? $opt{f} : 0.01);
my $WORKDIR = ($opt{w} ? $opt{w} : "/tmp/");
my $output = ($opt{o} ? $opt{o} : "");

# Check binaries
if( !(`$SNVMIX -h 2>&1` =~ m/^Syntax/) ) {
	die("ERROR: SNVMix2 binary not found, place in \$PATH or change \$SNVMIX = '$SNVMIX'\n");
}
if( !(`$ROCOPS -h 2>&1` =~ m/^Syntax/) ) {
	die("ERROR: ROCops script not found, place in \$PATH or change \$ROCOPS = '$ROCOPS'\n");
}

# get model, bQs and mQs values
my @bQs;
get_quality_vals($opt_bQs, \@bQs, "base");

my @mQs;
get_quality_vals($opt_mQs, \@mQs, "map");

my @models = split(/,/, $opt_model);

print STDERR "Evaluating SNVMix2 using all combinations of following parameters:\n";
print STDERR "\tBase qualities: ".join(" ", @bQs)."\n";
print STDERR "\tMap qualities: ".join(" ", @mQs)."\n";
print STDERR "\tModels:";
foreach my $model (@models) {
	if(exists($VALID_MODELS{$model})) { print STDERR " $model"; }
	else { die("\nERROR: model '$model' is not recognized as a valid model\n"); }
}
print STDERR "\n";


# Check train and test files
my %train_pos;
my %test_pos;
get_truth($train_set_file, \%train_pos);
get_truth($test_set_file, \%test_pos);


# Check there is no overlap between train and test sets
foreach my $pos (keys %train_pos) {
	if(exists($test_pos{$pos})) {
		die("ERROR: Position $pos is in both train and test sets\n");
	}
}

# Create test and train pileup files
if(! -w $WORKDIR) {
	die("ERROR: Work directory not writeable, \$WORKDIR = '$WORKDIR'\n");
}
$WORKDIR .= "/".$ENV{USERNAME}."-".int(rand(100000))."/";
if(!mkdir($WORKDIR)) {
	die("ERROR: Could not create working directory '$WORKDIR'\n");
}
print STDERR "Using $WORKDIR as work directory\n";

my $pileup_test = $WORKDIR."test.pileup";
my $pileup_train = $WORKDIR."train.pileup";
my $pileup_format = split_filter_pileups($pileup, $pileup_train, $pileup_test, \%train_pos, \%test_pos);


# Process respective models and get results for each run
my @results;
foreach my $model (@models) {
	if($VALID_MODELS{$model} == 1) {
		my $mQ = 0;
		foreach my $bQ (sort {$a <=> $b} @bQs) {
			print STDERR "Running: ${model}_mQ${mQ}_bQ${bQ}\n";
			my %run_result = get_snvmix_run_results($WORKDIR, $pileup_test, $pileup_train, \%test_pos, \%train_pos, $model, $bQ, $mQ, $pileup_format, $FDR);
			push @results, {%run_result};
		}
	} elsif($VALID_MODELS{$model} == 2) {
		my $bQ = 0;
		foreach my $mQ (sort {$a <=> $b} @mQs) {
			print STDERR "Running: ${model}_mQ${mQ}_bQ${bQ}\n";
			my %run_result = get_snvmix_run_results($WORKDIR, $pileup_test, $pileup_train, \%test_pos, \%train_pos, $model, $bQ, $mQ, $pileup_format, $FDR);
			push @results, {%run_result};
		}
	} elsif($VALID_MODELS{$model} == 3) {
		foreach my $mQ (sort {$a <=> $b} @mQs) {
			foreach my $bQ (sort {$a <=> $b} @bQs) {
				print STDERR "Running: ${model}_mQ${mQ}_bQ${bQ}\n";
				my %run_result = get_snvmix_run_results($WORKDIR, $pileup_test, $pileup_train, \%test_pos, \%train_pos, $model, $bQ, $mQ, $pileup_format, $FDR, "debug");
				push @results, {%run_result};
			}
		}
	} else {
		die("ERROR: Should not really be here, script corrupted?\n");
	}

}

# Calculate run(s) that gave max fmeasure
my $max_fmeas = -1;
my @max_fmeas_ids;
for(my $i = 0; $i < scalar(@results); $i++) {
	if($results[$i]{fmeas} > $max_fmeas) {
		$max_fmeas = $results[$i]{fmeas};
		$#max_fmeas_ids = -1;
	}
	if($results[$i]{fmeas} == $max_fmeas) {
		push @max_fmeas_ids, $i;
	}
}
if($max_fmeas == -1) {
	die("ERROR: unexpected error when finding model with max f-measure\n");
}
print STDERR "Max fmeasure = $max_fmeas with";
if(scalar(@max_fmeas_ids) == 1) {
	# max f-measure was not shared between runs
	print STDERR " model ".$results[$max_fmeas_ids[0]]{run_name}." THR = ".$results[$max_fmeas_ids[0]]{THR}."\n";
} else {
	# ah, several runs had same maximum f-measure
	print STDERR " models: \n";
	foreach my $run_id (@max_fmeas_ids) {
		print STDERR "\t".$results[$run_id]{run_name}." THR = ".$results[$run_id]{THR}."\n";
	}
}

if($output) {
	open(OUTPUT, ">$output") || die("ERROR: could not open output file $output for writing full scan results\n");
	print STDERR "Printing full scan output to $output\n";
	print OUTPUT join("\t","run_name","AUC","THR","pos_total","T","F","TP","FP","TN","FN","sens","spec","prec","fmeas")."\n";
	for(my $i = 0; $i < scalar(@results); $i++) {
	        print OUTPUT join("\t",$results[$i]{run_name},$results[$i]{AUC},$results[$i]{THR})."\t";
		print OUTPUT join("\t",$results[$i]{totalCalls},$results[$i]{T},$results[$i]{F})."\t";
		print OUTPUT join("\t",$results[$i]{TP},$results[$i]{FP},$results[$i]{TN},$results[$i]{FN})."\t";
		print OUTPUT join("\t",$results[$i]{sens},$results[$i]{spec},$results[$i]{prec},$results[$i]{fmeas})."\n";
	}
	close(OUTPUT);
}



###############
# SUBROUTINES #
###############
#my %run_result = get_snvmix_run_results($WORKDIR, $pileup_test, $pileup_train, \%test_pos, \%train_pos, "MB", 10, 20, 'm', 0.01);

sub get_snvmix_run_results() {
# Way this is called:
# %run_result = get_snvmix_run_results($WORKDIR, $pileup_test, $pileup_train, \%test_pos, \%train_pos, $model, $bQ, $mQ, $pileup_format, $FDR);
	my ($WORKDIR, $pileup_test, $pileup_train, $test_pos, $train_pos, $model, $bQ, $mQ, $pileup_format, $FDR) = @_;
	my $run_name = $model."_mQ".$mQ."_bQ".$bQ;

	# Sanity checks
	my $error = "";
	if(! -w $WORKDIR) { $error = "work directory is not writeable: '$WORKDIR'"; }
	if(! -r $pileup_test) { $error = "cannot read pileup_test file: '$pileup_test'"; }
	if(! -r $pileup_train) { $error = "cannot read pileup_train file: $pileup_train'"; }
	if(!($bQ =~ /^\d+$/)) { $error = "base quality is not a natural number: '$bQ'"; }
	if(!($mQ =~ /^\d+$/)) { $error = "map quality is not a natural number: '$mQ'"; }
	if(!($model =~ /^(mb|b|m|M|Mb|MB|SNVMix1)$/)) { $error = "model is not recognized: '$model'"; }
	if( !( $FDR <= 1 && $FDR >= 0) ) { $error = "FDR does not seem to be a decimal between 0 and 1: '$FDR'"; }
	if( $pileup_format ne 'm' && $pileup_format ne 's' ) { $error = "pileup format not recognized: '$pileup_format'"; }
	if($error) {
		print STDERR "ERROR trying to evaluate '$run_name': $error\n";
		return 0;
	}

	# Train SNVMix2
	my $param_file = $WORKDIR."/".$run_name.".model";
	`$SNVMIX -T -i $pileup_train -m $param_file -t $model -p $pileup_format -q $bQ -Q $mQ 2>/dev/null`;

	# Get calls for train set
	my $train_calls = `$SNVMIX -C -i $pileup_train -m $param_file -t $model -p $pileup_format -q $bQ -Q $mQ -f 2>/dev/null`;

	# Get threshold
	my $THR;
	my $AUC;
	my ($pipe_read, $pipe_write);
	my $pid = open2($pipe_read, $pipe_write, "$ROCOPS -f $FDR - /dev/null 2>&1");
	foreach my $train_call (split(/\n/, $train_calls)) {
		#chr1:878522	T	C	T:0,C:18,0.0000000000,0.0000035823,0.9999964177,3
		my ($pos, $ref, $nref, $data) = split(/\t/, $train_call);
		my ($ref_count, $nref_count, $pAA, $pAB, $pBB) = split(/,/, $data);
		if(exists($train_pos->{$pos})) {
			print $pipe_write ($pAB+$pBB)."\t".$train_pos->{$pos}."\n";
		} # else, pileup was not filtered
	}
	undef($train_calls);
	close($pipe_write);
	while(<$pipe_read>) {
		chomp;
		if(m/At FDR = $FDR, THR = (.+)$/) { $THR = $1; }
		if(m/AUC = (.+)$/) { $AUC = $1; }
	}
	close($pipe_read);
	if(!$THR) {
		print STDERR "ERROR trying to evaluate '$run_name': could not read THR from $ROCOPS script\n";
		return 0;
	}
	if(!$AUC) {
		print STDERR "ERROR trying to evaluate '$run_name': could not read AUC from $ROCOPS script\n";
		return 0;
	}

	# Get calls for test set
	my $test_calls = `$SNVMIX -C -i $pileup_test -m $param_file -t $model -p $pileup_format -q $bQ -Q $mQ -f 2>/dev/null`;

	# Get TP, TN, FP, FN, sens, prec, f-measure
	my $totalCalls = 0;
	my $T = 0; my $F = 0;
	my $TP = 0; my $FP = 0;
	my $TN = 0; my $FN = 0;
	foreach my $test_call (split(/\n/, $test_calls)) {
		#chr1:878522	T	C	T:0,C:18,0.0000000000,0.0000035823,0.9999964177,3
		my ($pos, $ref, $nref, $data) = split(/\t/, $test_call);
		my ($ref_count, $nref_count, $pAA, $pAB, $pBB) = split(/,/, $data);
		if(!exists($test_pos->{$pos})) {
			next;
		}
		my $truth = $test_pos->{$pos};
		my $call = $pAB+$pBB;
		if($truth) {
			$T++;
			if($call >= $THR) { $TP++; }
			else { $FN++; }
		} else {
			$F++;
			if($call >= $THR) { $FP++; }
			else { $TN++; }
		}
		$totalCalls++;
	}
	my $sens = "NaN";
	my $spec = "NaN";
	my $prec = "NaN";
	my $fmeas = "NaN";
	if( ($TP + $FN) ) { $sens = sprintf("%0.4f",$TP/($TP + $FN)); }
	if( ($TN + $FP) ) { $spec = sprintf("%0.4f",$TN/($TN + $FP)); }
	if( ($TP + $FP) ) { $prec = sprintf("%0.4f",$TP/($TP + $FP)); }
	if( ($TP + $FN) && ($TP + $FP) ) { $fmeas = sprintf("%0.4f",2*($prec*$sens)/($prec+$sens)); }

	#print join("\t",$run_name,$AUC,$THR,$totalCalls,$T,$F,$TP,$FP,$TN,$FN,$sens,$spec,$prec,$fmeas)."\n";

	return run_name => $run_name,
		AUC => $AUC,
		THR => $THR,
		totalCalls => $totalCalls,
		T => $T,
		F => $F,
		TP => $TP,
		FP => $FP,
		TN => $TN,
		FN => $FN,
		sens => $sens,
		spec => $spec,
		prec => $prec,
		fmeas => $fmeas;
}


sub get_quality_vals() {
	my $opt_Q = shift @_;
	my $Q = shift @_;
	my $type = shift @_;

	if($opt_Q =~ m/^(\d+)$/) {
		push @$Q, $1;
	} elsif($opt_Q =~ m/^(\d+):(\d+):(\d+)$/) {
		for(my $i = $1; $i <= $3; $i += $2) { push @$Q, $i; }
	} elsif($opt_Q =~ m/^(\d+),(\d+)/) {
		my @tmp = split(/,/, $opt_Q);
		my %check;
		foreach my $i (@tmp) {
			if(!($i =~ /^\d+$/)) {
				die("ERROR: invalid value in quality list: $i\n");
			}
			if(!exists($check{$i})) {
				push @$Q, $i;
				$check{$i}++;
			}
		}
	} else {
		die("ERROR: invalid format in $type quality list\n");
	}
}

sub split_filter_pileups() {
	my $pileup = shift @_;
	my $pileup_train = shift @_;
	my $pileup_test = shift @_;
	my $train_pos = shift @_;
	my $test_pos = shift @_;

	open(PILEUP_TEST,">$pileup_test") || die("ERROR: could not open test pileup file for writing ($pileup_test)\n");
	open(PILEUP_TRAIN,">$pileup_train") || die("ERROR: could not open test pileup file for writing ($pileup_train)\n");
	open(PILEUP,"<$pileup") || die("ERROR: could not open input pileup file for reading ($pileup)\n");
	my $line_number = 1;
	my $pileup_format = "";
	while(<PILEUP>) {
		chomp;
		s///;
		#if(m/^(\S+)\t(\d+)\t\w\t\d+\t@\S*\t@\S*\t@\S*$/) {
		if(m/^(\S+)\t(\d+)\t\w\t\d+\t@\S*\t@\S*\t@\S*/) {
			if(!$pileup_format) { $pileup_format = "m" }
			elsif($pileup_format ne 'm') { die("ERROR: line $line_number appears to change in format\n"); }
		} elsif(m/^(\S+)\t(\d+)\t\w\t\d+\t\S*\t\S*\t\S*$/) {
			if(!$pileup_format) { $pileup_format = "s" }
			elsif($pileup_format ne 's') { die("ERROR: line $line_number appears to change in format\n"); }
		} else {
			die("ERROR: line $line_number in pileup file $pileup fails format check\n");
		}
		my $pos = $1.":".$2;
		if(exists($train_pos->{$pos})) {
			print PILEUP_TRAIN "$_\n";
		} elsif(exists($test_pos->{$pos})) {
			print PILEUP_TEST "$_\n";
		}
		$line_number++;
	}
	close(PILEUP); close(PILEUP_TEST); close(PILEUP_TRAIN);
	return $pileup_format;
}

sub get_truth() {
	my $file = shift @_;
	my $pos = shift @_;

	open(INPUT, "<$file") || die("ERROR: Could not open 'truth' file for reading: $file\n");
	my $line_number = 1;
	while(<INPUT>) {
		chomp;
		s///;
		if(!m/^(\w+:\d+)\t([01])$/) {
			die("ERROR: file $file line $line_number does not comply with format: '$_'\n");
		}
		if(exists($pos->{$1})) {
			if($pos->{$1} != $2) {
				die("ERROR: file $file line $line_number provides a conflicting value for position $1\n");
			}
		} else {
			$pos->{$1} = $2;
		}
		$line_number++;
	}
	close(INPUT);
}



sub usage() {
        print "Syntax:\n$0\n";
        print << "EOF";
	-p <pileup file>
	-t <train set, chr:pos 0/1>
	-T <test set, chr:pos 0/1>
	-b <base Q, either fixed, comma separated list or start:inc:stop>
	-m <map Q, either fixed, comma separated list or start:inc:stop>
	-M <model code(s), e.g. MB,m,mb>
	-f <desired FDR>
	-o <full output file>

   	-h      this message
EOF
exit;
}



