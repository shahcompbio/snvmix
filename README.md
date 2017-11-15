*** SNVMix2 ***

This version of SNVMix2 has has been tested under Linux and Mac OS X.

To build:

> unzip -x ../SNVMix2-v{VERSION}.zip
> cd SNVmix2-v{VERSION}/
> make

Binary file is called "SNVMix2", for help run with -h flag

> ./SNVMix2 -h

You can copy that file to a preferred location.

SNVMix2 defaults to reading from standard input if flag '-i' is not specified,
so you can use it in pipe mode next to a MAQ or samtools pileup command, saving
storage space. Same applies for standard output and the "-o" flag while Classifying.

Parameter '-m <file>' is used to read the model parameters when Classifying (-C)
or write them in Train mode (-T)

Different models for base and mapping qualities are available using the "-t" flag,
SNVMix1 mode can be accessed by selecting "-t SNVMix1"

Pileup file should be generated with base and mapping qualities but without consensus,
such as the one obtained when running, with samtools-0.1.4:
	samtools pileup -s -l <in_list.txt> -f <ref.fa> <in.bam>

or with maq-0.7.0:
	maq pileup -v <ref.bfa> <in.map>

SNVMix will output 4 columns:
1: coordinate in "seq:pos" format
2: reference base
3: non-reference base
4: comma separated field:
			REF:#, NREF:#, p(AA), p(AB), p(BB), maxP

	REF:#	reference base and number of occurrences that passed quality settings
	NREF:#	non-reference base and number of occurrences that passed quality settings
	p(AA)	probability assigned to homozygous to reference
	p(AB)	probability assigned to heterozygous genotype
	p(BB)	probability assigned to homozygous to the non-reference
	maxP 	class with max probability (1=AA, 2=AB, 3=BB)


A perl script is provided to filter SNVMix2 result according to a specified
threshold, h
elp on this script can be found running:

> ./misc/snvmix2summary.pl -h

This script can filter candidate SNVs using two methods specified by the '-c <TYPE>'
flag:

'-c 2'	Will consider only two classes, either homozygote for the reference allele
	( p(AA) ) or not ( p(AB) + p(BB) )

'-c 3'	Will consider the three clases p(AA), p(AB), p(BB) separately


The output still still retains the unmodified probability values, so p_AB can still
be distinguished from p_BB in case '-c 2' is used.





Notes:

A working gcc compiler is needed and under Linux libc >= 4.6.27 is required.

Be careful if you copy & paste pileup files, as some text editors tend to change
"tabs" for spaces, which is needed as the field delimiter. This can cause problems
in the parser.

//
Comments and Questions:
Sohrab Shah	sshah@bccrc.ca
Rodrigo Goya	rgoya@bcgsc.ca

//
LICENSE: this software is distributed under the MIT license
