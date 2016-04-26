/* The MIT License

   Copyright (c) 2009, by Sohrab Shah <sshah@bccrc.ca> and Rodrigo Goya <rgoya@bcgsc.ca>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/*
C Implementation of SNVMix2
*/

#define VERSION "0.11.8-r4"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SNVMix2.h"

#define START_QNUM 1000

int main(int argc, char **argv) {
	
	param_struct params;// = {NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0};

	initSNVMix(argc, argv, &params);

	if(params.classify || params.filter) {
		if(params.input_type == Q_PILEUP)
			snvmixClassify_qualities(&params);
		else if(params.input_type == M_PILEUP || params.input_type == S_PILEUP)
			snvmixClassify_pileup(&params);
	} else if(params.train) {
		if(params.input_type == M_PILEUP || params.input_type == S_PILEUP)
			snvmixTrain_pileup(&params);
		else if(params.input_type == Q_PILEUP) {
			fprintf(stderr,"Sorry, Training with quality-type input is not yet implemented\n");
			exit(1);
		}
	}

}

/*
	void snvmixClassify_qualities
	Arguments:
		*params : pointer to model parameters and command line options
	Function:
		Reads a pilep-style file that has been preprocessed to include the quality of
		calls being the reference allele and the maping quality both in decimal values.

		For each line it evaluates de model according to the parameters in *params and
		outputs the corresponding SNV prediction values.

		This function may come in handy when filtering of calls is done by a
		different program or the base-calling and maping quality information is not
		in a pileup/phred format.
	
		The file read is a tab-separated file where the columns are:
			-  target sequence (e.g. chromosome)
			-  sequence position
			-  reference base
			-  non-reference base
			-  comma separated probabilities of the call being the reference
			-  comma separated mapping qualities, one for each base call

*/
void snvmixClassify_qualities(param_struct *params) {
	FILE *fptrIN, *fptrOUT;

	char *line = NULL;
	size_t line_len = 0, prev_len = 0;
	int read_len = 0;

	int qual_num = 0, col_num = 0;

	char *col, *qual;
	char *col_sep = "\t", *col_str, *col_ptr;
	char *qual_sep = ",", *qual_str, *qual_ptr;
	char *chr_pos, *ref, *nref;	

	long double *bQ, *mQ, *tmpQ;
	int depth = 0, maxp;
	long int depth_allocated = 0;

	long double bsum = 0, msum = 0, b_AA = 0, b_AB = 0, b_BB = 0, mu[3], notmu[3], pi[3], log_pi[3], z = 0;
	int i;

	for(i = 0; i < 3; i++) {
		mu[i] = params->mu[i];
		notmu[i] = 1 - mu[i];
		pi[i] = params->pi[i];
		log_pi[i] = logl(pi[i]);
	}
	fptrIN = params->input;
	fptrOUT = params->output;


	if( (bQ = malloc(START_QNUM * sizeof (long double))) == NULL) {
		perror("ERROR allocating initial space for bQ"); exit(1);
	}
	if( (mQ = malloc(START_QNUM * sizeof (long double))) == NULL) {
		perror("ERROR allocating initial space for mQ"); exit(1);
	}
	depth_allocated = START_QNUM;
#if defined __linux__ || defined _GNU_SOURCE
	while( (read_len = getline(&line,&line_len,fptrIN)) != -1 ) {
#endif
#ifdef __APPLE__
	while( (line = fgetln(fptrIN,&line_len)) != NULL ) {
		line[line_len-1] = '\0';
#endif
		depth = 0;
		col_ptr = NULL;
		for(col_num = 0, col_str = line; ; col_num++, col_str = NULL) {
			col = strtok_r(col_str, "\t", &col_ptr);
			if(col == NULL) {
				break;
			}
			if(col_num == 0) {
				chr_pos = col;
			} else if(col_num == 1) {
				ref = col;
			} else if(col_num == 2) {
				nref = col;
			} else if(col_num == 3) {
				for(qual_num = 0, qual_str = col; ; qual_num++, qual_str = NULL) {
					qual = strtok_r(qual_str, ",", &qual_ptr);
					if(qual == NULL) {
						break;
					}
					if(qual_num >= depth_allocated) {
						if( (bQ = realloc(bQ, sizeof(long double) * (depth_allocated + START_QNUM))) == NULL) {
							perror("ERROR allocating bQ"); exit(1);
						}
						if( (mQ = realloc(mQ, sizeof(long double) * (depth_allocated + START_QNUM))) == NULL) {
							perror("ERROR allocating mQ"); exit(1);
						}
						depth_allocated = depth_allocated + START_QNUM;
					}
					bQ[qual_num] = atof(qual);
				}
				depth = qual_num;
			} else if(col_num == 4) {
				for(qual_num = 0, qual_str = col; ; qual_num++, qual_str = NULL) {
					qual = strtok_r(qual_str, ",", &qual_ptr);
					if(qual == NULL) {
						break;
					}
					if(qual_num >= depth_allocated) {
						fprintf(stderr, "FATAL ERROR: should not have more mapping qualities than base qualities\n");
						exit(1);
					}
					mQ[qual_num] = atof(qual);
				}
				if(depth != qual_num) {
					fprintf(stderr, "FATAL ERROR: should not have more base qualities than mapping qualities\n");
					exit(1);
				}
			}
		}


		b_AA = 0;
		b_AB = 0;
		b_BB = 0;
    		//b = sum(log(notm*0.5 + m.*(yr.*mur+notyr.*notmur)),2);
		for(qual_num = 0; qual_num < depth; qual_num++) {
			b_AA = b_AA + logl((1-mQ[qual_num])*0.5 + mQ[qual_num]*(bQ[qual_num]*mu[0]+(1-bQ[qual_num])*notmu[0]));
			b_AB = b_AB + logl((1-mQ[qual_num])*0.5 + mQ[qual_num]*(bQ[qual_num]*mu[1]+(1-bQ[qual_num])*notmu[1]));
			b_BB = b_BB + logl((1-mQ[qual_num])*0.5 + mQ[qual_num]*(bQ[qual_num]*mu[2]+(1-bQ[qual_num])*notmu[2]));
		}
		b_AA = expl(b_AA + log_pi[0]);
		b_AB = expl(b_AB + log_pi[1]);
		b_BB = expl(b_BB + log_pi[2]);
		z = b_AA + b_AB + b_BB;
		if(!z) { z = 1; }
		b_AA = b_AA / z;
		b_AB = b_AB / z;
		b_BB = b_BB / z;
		maxp = 1;
		z = b_AA;
		if( b_AB > z) { maxp = 2; z = b_AB; }
		if( b_BB > z) { maxp = 3; }

		fprintf(fptrOUT,"%s\t%s\t%s",chr_pos, ref, nref);
		fprintf(fptrOUT,"\t%0.10Lf,%0.10Lf,%0.10Lf,%d\n",b_AA,b_AB,b_BB,maxp);
	}
	fclose(fptrIN);
	fclose(fptrOUT);
	free(line);
}


/*
	void snvmixClassify_pileup
	Arguments:
		*params : pointer to model parameters and command line options
	Function:
		Reads a MAQ or Samtools pileup format. For required formatting we refer the user
		to the corresponding websites:
			http://maq.sourceforge.net/
			http://samtools.sourceforge.net/

		In general the format of both files can be described as a tab separated file
		where columns are:
			-  target sequence (e.g. chromosome)
			-  sequence position
			-  reference base
			-  depth
			-  stream of base calls
			-  stream of base call PHRED-scaled qualities
			-  stream of mapping PHRED-scaled qualities

		Note: when handling pileup files, care should be taken when doing copy-paste
		operations not to transform 'tabs' into streams of spaces.

		For each line it evaluates de model according to the parameters in *params and
		outputs the corresponding SNV prediction values.

		Predictions can be output either only for positions where when non-reference values
		are observed or, when run with '-f' flag, for all positions. The -f flag is useful
		when the resulting calls are being compared or joined accros different datasets
		and all predictions need to be present.
*/
void snvmixClassify_pileup(param_struct *params) {
//MAQ:
//	1       554484  C       1752    @,,.,.,, @xxx @xxx xxxxx
//SAM:
//	1       554484  C       1752    ,,.,.,,	xxx xxxx

	FILE *fptrIN, *fptrOUT;

	char *line = NULL;
	size_t line_len = 0, prev_len = 0;
	int read_len = 0;
	int col_num = 0;
	long int line_num = 0;

	char *col, *qual;
	char *col_sep = "\t", *col_str, *col_ptr;
	char *qual_sep = ",", *qual_str, *qual_ptr;
	char *chr, *pos, *ref, nref, call, nref_skip = 'N';

	char *bQ, *mQ;
	int *calls, *tmpQ;
	int depth = 0, qual_num,  maxp;
	int depth_allocated = 0;

	long double bsum = 0, msum = 0, b_AA = 0, b_AB = 0, b_BB = 0, mu[3], notmu[3], pi[3], log_pi[3], z = 0;
	int nref_count = 0, ref_count = 0;
	long double Y, not_Y;
	long double phred[PHRED_MAX + 1];
	int i, call_i;

	//initPhred(phred, PHRED_MAX+1);
	initPhred(phred, PHRED_MAX);

	for(i = 0; i < 3; i++) {
		mu[i] = params->mu[i];
		notmu[i] = 1 - mu[i];
		pi[i] = params->pi[i];
		log_pi[i] = logl(pi[i]);
	}
	fptrIN = params->input;
	fptrOUT = params->output;
	if(params->full)
		nref_skip = -2;


	if( (calls = malloc(START_QNUM * sizeof (int))) == NULL) {
		perror("ERROR allocating initial space for calls"); exit(1);
	}
	depth_allocated = START_QNUM;
#if defined __linux__ || defined _GNU_SOURCE
	while( (read_len = getline(&line,&line_len,fptrIN)) != -1 ) {
#endif
#ifdef __APPLE__
	while( (line = fgetln(fptrIN,&line_len)) != NULL ) {
		line[line_len-1] = '\0';
#endif
		line_num++;
		depth = 0;
		nref = 0;
		col_ptr = NULL;
		for(col_num = 0, col_str = line; nref != nref_skip ; col_num++, col_str = NULL) {
			col = strtok_r(col_str, "\t", &col_ptr);
			if(col == NULL) {
				break;
			}
			if(col_num == 0) {
				chr = col;
			} else if(col_num == 1) {
				pos = col;
			} else if(col_num == 2) {
				ref = col;
			} else if(col_num == 3) {
				depth = atoi(col);
				if(depth > depth_allocated) {
					if( (calls = realloc(calls, sizeof (int) * (depth + START_QNUM))) == NULL) {
						perror("ERROR allocating space for calls"); exit(1);
					}
					depth_allocated = depth + START_QNUM;
				}
			} else if(col_num == 4) {
				if(params->input_type == M_PILEUP)
					col = col+sizeof(char);
					snvmixGetCallString(col, calls, depth, &nref);
			} else if(col_num == 5) {
				bQ = col;
				// If it's a MAQ pileup, we need to skip the @ sign
				if(params->input_type == M_PILEUP)
					bQ = bQ + sizeof(char);
				for(call_i = 0; call_i < depth; call_i++)
					bQ[call_i] = bQ[call_i]-33;
			} else if(col_num == 6) {
				mQ = col;
				// If it's a MAQ pileup, we need to skip the @ sign
				if(params->input_type == M_PILEUP)
					mQ = mQ + sizeof(char);
				for(call_i = 0; call_i < depth; call_i++)
					mQ[call_i] = mQ[call_i]-33;
			}
		}
		if(col_num != 7) {
			fprintf(stderr, "ERROR: Pileup line %d has wrong number of columns (has %d expected 7, was flag '-s' used in samtools to create pileup?)\n", line_num, col_num);
			exit(1);
		}
		// If we quit the for because no nref was found, skip this too, next line
	nref = snvmixFilterCalls(calls,depth,bQ,mQ,params);
	nref_count = 0;
	ref_count = 0;
	for(qual_num = 0; qual_num < depth; qual_num++) {
		//if( snvmixSkipCall(calls,qual_num,params,bQ,mQ) )
		if(calls[qual_num] == -1)
			continue;
		if(calls[qual_num] == 1)
			ref_count++;
		if(calls[qual_num] == nref)
			nref_count++;
	}
if(params->filter) {
	fprintf(fptrOUT,"%s:%s\t%s:%d\t%c:%d\t", chr, pos, ref, ref_count, nref, nref_count);
	for(qual_num = 0; qual_num < ref_count; qual_num++)
		fprintf(fptrOUT,",");
	for(qual_num = 0; qual_num < nref_count; qual_num++)
		fprintf(fptrOUT,"%c",nref);
	fprintf(fptrOUT,"\n");
} else {
		if(nref == nref_skip)
			continue;
		//nref = snvmixFilterCalls(calls,depth,bQ,mQ,params);
		//if(nref == nref_skip)
		//	continue;
		b_AA = 0;
		b_AB = 0;
		b_BB = 0;
		nref_count = 0;
		for(qual_num = 0; qual_num < depth; qual_num++) {
			//if( snvmixSkipCall(calls,qual_num,params,bQ,mQ) )
			if(calls[qual_num] == -1)
				continue;
		// from matlab:
		// b = sum(log(notm*0.5 + m.*(yr.*mur+notyr.*notmur)),2);
			if(calls[qual_num] == 1) {
				// If REF then
				Y = phred[(unsigned char) bQ[qual_num]];
				not_Y = 1-phred[(unsigned char) bQ[qual_num]];
			} else {
				nref_count++;
				// If NREF then
				Y = (1-phred[(unsigned char) bQ[qual_num]])/3;
				not_Y = phred[(unsigned char) bQ[qual_num]] + 2*(1-phred[(unsigned char)bQ[qual_num]])/3;
			}
			b_AA = b_AA + logl((1-phred[(unsigned char) mQ[qual_num]])*0.5 + phred[(unsigned char) mQ[qual_num]] * (Y * mu[0] + not_Y * notmu[0]));
			b_AB = b_AB + logl((1-phred[(unsigned char) mQ[qual_num]])*0.5 + phred[(unsigned char) mQ[qual_num]] * (Y * mu[1] + not_Y * notmu[1]));
			b_BB = b_BB + logl((1-phred[(unsigned char) mQ[qual_num]])*0.5 + phred[(unsigned char) mQ[qual_num]] * (Y * mu[2] + not_Y * notmu[2]));
		}
		// Check if any non-references actually passed the filter
		//if(!nref_count)
		//	continue;
		b_AA = expl(b_AA + log_pi[0]);
		b_AB = expl(b_AB + log_pi[1]);
		b_BB = expl(b_BB + log_pi[2]);
		z = b_AA + b_AB + b_BB;
		if(!z) z = 1;
		b_AA = b_AA / z;
		b_AB = b_AB / z;
		b_BB = b_BB / z;
		maxp = 1;
		z = b_AA;
		if( b_AB > z) { maxp = 2; z = b_AB; }
		if( b_BB > z) { maxp = 3; }


		fprintf(fptrOUT,"%s:%s\t%s\t%c\t%s:%d,%c:%d,",chr,pos, ref, nref, ref, ref_count,  nref, nref_count);
		fprintf(fptrOUT,"%0.10Lf,%0.10Lf,%0.10Lf,%d\n",b_AA,b_AB,b_BB,maxp);
}
	}
	fclose(fptrIN);
	fclose(fptrOUT);
	free(line);
}

/*
	void snvmixGetCallString
	Arguments:
		*col : pointer to the current file-line in memory
		*calls : pointer to array where we will store the final base-calls
		depth : number of base reads for this position
		*nref : the observed NREF value will be placed here (-1 if none was found)

	Function:
		This will parse column 5 of the pileup file, which contains the
		base calls and will fill the array "calls" with:
			1 : when a reference call was made
		       -1 : when an invalid value is seen ('N', 'n', '*') 
		   [ACTG] : when a non-reference call was made

		

*/
void snvmixGetCallString(char *col, int *calls, int depth, char *nref) {
	int i;
	int call_i = 0, call_skip_num = 0;
	char call_skip_char[10];
	for(i = 0 ; call_i < depth; i++) {
		switch(col[i]){
		case ',':
		case '.':
			calls[call_i] = 1;
			call_i++;
			break;
		case 'a': case 'A':
		case 't': case 'T':
		case 'g': case 'G':
		case 'c': case 'C':
			//calls[call_i] = 0;
			calls[call_i] = toupper(col[i]);
			call_i++;
			//if(*nref == 0)
				//*nref = toupper(*(col+sizeof(char)*i));
			if(*nref == 0)
				*nref = toupper(col[i]);
				break;
		case 'N':
		case 'n':
		case '*':
		case '<': // These values are inserted by samtools > 0.1.8
		case '>': // to represent "skipping" events
			// Not comparable values, we ignore them, but need to
			// keep count to compare against the number of mapping qualities
			calls[call_i] = -1;
			call_i++;
			break;
		case '$':
			// End of a read, just ignore it
			break;
		case '^':
			// Start of read, ignore it and skip the next value (not base-related)
			i++;
			break;
		case '+':
		case '-':
			// Eliminate:
			//		+2AA
			//		-3AAA
			// Start skiping the sign
			i++;
			for(call_skip_num = 0; col[i] <= '9' && col[i] >= '0'; call_skip_num++, i++) {
				//call_skip_char[call_skip_num] = call;
				call_skip_char[call_skip_num] = col[i];
				//i++;
			}
			// we have the number string in call_skip_char, lets terminate it
			call_skip_char[call_skip_num] = '\0';
			// i has been updated to first letter, just need to jump the rest of them
			i = i+atoi(call_skip_char)-1;
			break;
		default:
			fprintf(stderr,"ERROR: problem reading pileup file calls (%c)\n",col[i]);
			exit(1);
			break;
		}
	}
	// If no nref was found, don't bother parsing the other columns, make the for fail with -1
	if(!*nref)
		*nref = -1;
}

/*
	int snvmixFilterCalls
	Arguments:
		*calls : array built by snvmixGetCallString containing
			1 : when a reference call was made
		       -1 : when an invalid value is seen ('N', 'n', '*') 
		   [ACTG] : when a non-reference call was made
		depth : number of calls for the current position
		*bQ : base qualities observed for each call
		*mQ : map quality for each call
		*params : parameter structure, get thresholding data
	Function:
		For each valid call read in the position this function will apply
		thresholding according to the type selected (-t flag) and the bQ (-q)
		and mQ (-Q) thresholds provided.

		Any base-call that does not pass thresholds will be switched from it's
		current value in *calls to -1;

		The most prevalent NREF after filtering will be determined and returned.
	
*/
int snvmixFilterCalls(int *calls, int depth, char *bQ, char *mQ, param_struct *params) {
	int qual_num, nref_counts[5];
	nref_counts[0] = 0;
	nref_counts[1] = 0;
	nref_counts[2] = 0;
	nref_counts[3] = 0;
	nref_counts[4] = 0;
	int max_nref = N;

	char nref = 0;

	for(qual_num = 0; qual_num < depth; qual_num++) {
		if( snvmixSkipCall(calls,qual_num,params,bQ,mQ) ) {
			calls[qual_num] = -1;
		} else {
			nref_counts[0]++;
			switch(calls[qual_num]) {
				case 'A':
					nref_counts[A]++;
					break;
				case 'C':
					nref_counts[C]++;
					break;
				case 'G':
					nref_counts[G]++;
					break;
				case 'T':
					nref_counts[T]++;
					break;
				case 1:
				case -1:
					nref_counts[0]--;
					break;
				default:
					fprintf(stderr,"ERROR: unknown call base\n");
					exit(1);
			}
		}
	}
	if(nref_counts[0]) {
		max_nref = A;
		if(nref_counts[max_nref] < nref_counts[C])
			max_nref = C;
		if(nref_counts[max_nref] < nref_counts[G])
			max_nref = G;
		if(nref_counts[max_nref] < nref_counts[T])
			max_nref = T;
	} else {
		//return -1;
	}
	for(qual_num = 0; qual_num < depth; qual_num++) {
		if(calls[qual_num] == 1)
			continue;
		if(calls[qual_num] != base_code[max_nref])
			calls[qual_num] = -1;
	}
	return base_code[max_nref];
}

/*
	int snvmixSkipCall
	Arguments:
		*calls : array built by snvmixGetCallString containing
			1 : when a reference call was made
		       -1 : when an invalid value is seen ('N', 'n', '*') 
		   [ACTG] : when a non-reference call was made
		qual_num : call number that is being evaluated
		*params : parameter structure, get thresholding data
		*bQ : base qualities observed for each call
		*mQ : map quality for each call
	Function:
		Evalates quality values in each of the filtering models

		Returns 1 if the calls[qual_num] needs to be filtered out
		Returns 0 otherwise
*/
int snvmixSkipCall(int *calls, int qual_num, param_struct *params, char *bQ, char *mQ) {
	if(calls[qual_num] == -1)
		return(1);
	if(params->filter_type == TYPE_mb) {
		if(bQ[qual_num] <= params->bQ || mQ[qual_num] <= params->mQ)
			return(1);
	} else if(params->filter_type == TYPE_b) {
		if(bQ[qual_num] <= params->bQ)
			return(1);
	} else {
		if(mQ[qual_num] <= params->mQ)
			return(1);
		if(params->filter_type == TYPE_m) {
			// Use mapping as surrogate
			bQ[qual_num] = mQ[qual_num];
			// Make mapping one
			mQ[qual_num] = (char) PHRED_MAX;
		} else if(params->filter_type == TYPE_M) {
			// Mapping passed, make it one
			mQ[qual_num] = (char) PHRED_MAX;
		} else if(params->filter_type == TYPE_Mb) {
			// Nothing special here
		} else if(params->filter_type == TYPE_MB || params->filter_type == TYPE_SNVMix1) {
			if(bQ[qual_num] <= params->bQ)
				return(1);
			if(params->filter_type == TYPE_SNVMix1) {
				bQ[qual_num] = (char) PHRED_MAX;
				mQ[qual_num] = (char) PHRED_MAX;
		}	}
	}
	return(0);
}



void resetParams(param_struct *params) {
	params->input = stdin;
	params->output = stdout;
	params->inputfile = NULL;
	params->outputfile = NULL;
	params->modelfile = NULL;
	params->filter_type = 0;
	params->train = 0;
	params->classify = 0;
	params->filter = 0;
	params->full = 0;
	params->input_type = S_PILEUP; // 0 = processed, 1 = maq pileup, 2 = sam pileup
	params->mu[3];
	params->pi[3];
	params->max_iter = 1000;
	params->bQ = 19;
	params->mQ = 19;
	params->debug = 0;
	params->trainP.alpha[0] = 1000;
	params->trainP.alpha[1] = 500;
	params->trainP.alpha[2] = 100;
	params->trainP.beta[0] = 100;
	params->trainP.beta[1] = 500;
	params->trainP.beta[2] = 1000;
	params->trainP.delta[0] = 10;
	params->trainP.delta[1] = 1;
	params->trainP.delta[2] = 1;
	params->trainP.param_file = NULL;
	params->trainP.max_iter = 100;
	params->trainP.bQ = NULL;
	params->trainP.mQ = NULL;
	params->trainP.calls = NULL;
}

void initSNVMix(int argc, char **argv, param_struct *params) {
	char c;
	resetParams(params);
	while ((c = getopt (argc, argv, "hDTCFfi:o:m:t:p:q:Q:a:b:d:M:")) != -1) {
		switch (c)
		{
			case 'm':
				params->modelfile = optarg;
				break;
			case 'i':
				params->inputfile = optarg;
				break;
			case 'o':
				params->outputfile = optarg;
				break;
			case 'T':
				params->train = 1;
				break;
			case 'C':
				params->classify = 1;
				break;
			case 'F':
				params->filter = 1;
				break;
			case 'f':
				params->full = 1;
				break;
			case 'p':
				if(*optarg == 'q') {
					params->input_type = Q_PILEUP;
				} else if(*optarg == 'm') {
					params->input_type = M_PILEUP;
				} else if(*optarg == 's') {
					params->input_type = S_PILEUP;
				} else {
					fprintf(stderr,"ERROR: unknown pileup format %c\n",*optarg);
					exit(1);
				}
				break;
			case 't':
				if(strcmp("mb",optarg) == 0)
					params->filter_type = TYPE_mb;
				else if(strcmp("m",optarg) == 0)
					params->filter_type = TYPE_m;
				else if(strcmp("b",optarg) == 0)
					params->filter_type = TYPE_b;
				else if(strcmp("M",optarg) == 0)
					params->filter_type = TYPE_M;
				else if(strcmp("Mb",optarg) == 0)
					params->filter_type = TYPE_Mb;
				else if(strcmp("MB",optarg) == 0)
					params->filter_type = TYPE_MB;
				else if(strcmp("SNVMix1",optarg) == 0)
					params->filter_type = TYPE_SNVMix1;
				else {
					fprintf(stderr,"ERROR: filter type '%s' not recognized\n",optarg);
					exit(1);
				}
				break;
			case 'q':
				params->bQ = atoi(optarg);
				if( params->bQ < 0 || params->bQ >= PHRED_MAX )  {
					fprintf(stderr,"ERROR: quality threshold value Q%d out of range\n",params->bQ);
					exit(1);
				}
				break;
			case 'Q':
				params->mQ = atoi(optarg);
				if( params->mQ < 0 || params->mQ >= PHRED_MAX )  {
					fprintf(stderr,"ERROR: mapping threshold value Q%d out of range\n",params->mQ);
					exit(1);
				}
				break;
			case 'h':
				usage(argv[0]);
				break;
			case 'D':
				params->debug = 1;
				break;
			case 'a': /* Alpha parameters */
				if(sscanf(optarg, "%Lf,%Lf,%Lf", &(params->trainP.alpha[0]), &(params->trainP.alpha[1]), &(params->trainP.alpha[2])) != 3) {
					fprintf(stderr, "ERROR: could not read alpha parameters, expecting: #,#,#\n");
					exit(1);
				}
				break;
			case 'b': /* Beta parameters */
				if(sscanf(optarg, "%Lf,%Lf,%Lf", &(params->trainP.beta[0]), &(params->trainP.beta[1]), &(params->trainP.beta[2])) != 3) {
					fprintf(stderr, "ERROR: could not read beta parameters, expecting: #,#,#\n");
					exit(1);
				}
				break;
			case 'd': /* Delta parameters */
				if(sscanf(optarg, "%Lf,%Lf,%Lf", &(params->trainP.delta[0]), &(params->trainP.delta[1]), &(params->trainP.delta[2])) != 3) {
					fprintf(stderr, "ERROR: could not read delta parameters, expecting: #,#,#\n");
					exit(1);
				}
				break;
			case 'M':
				params->trainP.param_file = optarg;
				break;
			default:
				fprintf(stderr,"Unknown option\n");
				usage(argv[0]);
				break;
			}
	}
	// Decide if we're going to train or classify
	if(params->filter) {
		params->train = 0;
		params->classify = 0;
	}
	if(params->train && params->classify) {
		fprintf(stderr,"ERROR: choose either train or classify\n");
		exit(1);
	} else if(!params->train && !params->classify && !params->filter) {
		fprintf(stderr,"Train/Classify/Filter not selected, picking default: Classify\n");
		params->classify = 1;
	}
	if( params->train ) {
		if( params->trainP.param_file ) {
			readTrainingParams(params);
		}
	}
	if( params->modelfile != NULL ) {
		if( params->classify ) {
			readModel(params);
		}
	} else if(params->classify) {
		//fprintf(stderr, "ERROR: need to specify a model file\n");
		//exit(1);
		params->mu[0] = 0.995905287891696078261816182930;
		params->mu[1] = 0.499569401000925172873223800707;
		params->mu[2] = 0.001002216846753116409260431219;
		params->pi[0] = 0.672519580755555956841362785781;
		params->pi[1] = 0.139342327124010650907237618412;
		params->pi[2] = 0.188138092120433392251399595807;
		fprintf(stderr,"Parameter file not given, using for mQ35 bQ10:\n");
		fprintf(stderr,"\tmu[0] = %Lf\n", params->mu[0]);
		fprintf(stderr,"\tmu[1] = %Lf\n", params->mu[1]);
		fprintf(stderr,"\tmu[2] = %Lf\n", params->mu[2]);
		fprintf(stderr,"\tpi[0] = %Lf\n", params->pi[0]);
		fprintf(stderr,"\tpi[1] = %Lf\n", params->pi[1]);
		fprintf(stderr,"\tpi[2] = %Lf\n", params->pi[2]);
	}
	if( params->inputfile != NULL) {
		params->input = fopen(params->inputfile, "r");
		if(params->input == NULL) {
			perror("ERROR: could not open input file");
			exit(1);
		}
	} else {
		params->input = stdin;
	}
	if( params->outputfile != NULL ) {
		params->output = fopen(params->outputfile, "w");
		if(params->output == NULL) {
			perror("ERROR: could not open output file");
			exit(1);
		}
	} else {
		params->output = stdout;
	}
}

void readTrainingParams(param_struct *params) {
	FILE *fptrPARAM;
	char *line = NULL;
	size_t line_len = 0, read_len = 0;
	char *param, *param_ptr;
	int i;

	fptrPARAM=fopen(params->trainP.param_file, "r");
	if(fptrPARAM == NULL) {
		perror("ERROR: could not open training parameter file");
		exit(1);
	}

#if defined __linux__ || defined _GNU_SOURCE
	while( (read_len = getline(&line,&line_len,fptrPARAM)) != -1 ) {
#endif
#ifdef __APPLE__
	while( (line = fgetln(fptrPARAM,&line_len)) != NULL ) {
		line[line_len-1] = '\0';
#endif
		param = strtok_r(line, " ", &param_ptr);
		// All parameters come in triplets
		if(strcmp("alpha", param) == 0) {
			for(i = 0; i < 3; i++) {
				param = strtok_r(NULL, " ", &param_ptr);
				if(param == NULL) {
					fprintf(stderr,"ERROR: missing value #%d for alpha\n", i);
					exit(1);
				}
				params->trainP.alpha[i] = atof(param);
			}
		} else if(strcmp("beta", param) == 0) {
			for(i = 0; i < 3; i++) {
				param = strtok_r(NULL, " ", &param_ptr);
				if(param == NULL) {
					fprintf(stderr,"ERROR: missing value #%d for beta\n", i);
					exit(1);
				}
				params->trainP.beta[i] = atof(param);
			}
		} else if(strcmp("delta", param) == 0) {
			for(i = 0; i < 3; i++) {
				param = strtok_r(NULL, " ", &param_ptr);
				if(param == NULL) {
					fprintf(stderr,"ERROR: missing value #%d for delta\n", i);
					exit(1);
				}
				params->trainP.delta[i] = atof(param);
			}
		} else {
			fprintf(stderr,"Ignoring unknown parameter:\n",line);
		}
	}
	fclose(fptrPARAM);
	free(line);
}

void readModel(param_struct *params) {
	FILE *fptrPARAM;
	char *line = NULL;
	size_t line_len = 0, read_len = 0;
	char *param, *param_ptr;
	int i;

	fptrPARAM=fopen(params->modelfile, "r");
	if(fptrPARAM == NULL) {
		perror("ERROR: could not open parameter file");
		exit(1);
	}

#if defined __linux__ || defined _GNU_SOURCE
	while( (read_len = getline(&line,&line_len,fptrPARAM)) != -1 ) {
#endif
#ifdef __APPLE__
	while( (line = fgetln(fptrPARAM,&line_len)) != NULL ) {
		line[line_len-1] = '\0';
#endif
		param = strtok_r(line, " ", &param_ptr);
		if(strcmp("mu", param) == 0) {
			fprintf(stderr,"Setting values for mu:\n");
			// Need to read three values
			for(i = 0; i < 3; i++) {
				param = strtok_r(NULL, " ", &param_ptr);
				if(param == NULL) {
					fprintf(stderr,"ERROR: missing value #%d for mu\n", i);
					exit(1);
				}
				params->mu[i] = atof(param);
				fprintf(stderr,"\tmu[%d] = %Lf\n",i, params->mu[i]);
			}
		} else if(strcmp(param, "pi") == 0) {
			fprintf(stderr,"Setting valus for pi:\n");
			for(i = 0; i < 3; i++) {
				param = strtok_r(NULL, " ", &param_ptr);
				if(param == NULL) {
					fprintf(stderr,"ERROR: missing value #%d for pi\n", i);
					exit(1);
				}
				params->pi[i] = atof(param);
				fprintf(stderr,"\tpi[%d] = %Lf\n",i, params->pi[i]);
			}
		} else {
			fprintf(stderr,"Ignoring unknown parameter:\n",line);
		}
	}
	fclose(fptrPARAM);
	free(line);
}

void initPhred(long double *phredTable, int elem) {
	int i;
	for(i = 0; i < elem; i++) {
		phredTable[i] = 1-powl(10,(-(long double)i/10));
}
	phredTable[i] = 1;
}


void usage(char *selfname) {
	param_struct default_params;
	resetParams(&default_params);
	fprintf(stderr,"Syntax:\n\t%s\t-m <modelfile> [-i <infile>] [-o <outfile>] [-T | -C | -F] [-p < q | m | s >] [-t < mb | m | b | M | Mb | MB | SNVMix1>] [-q <#>] [-Q <#>] [-a/-b/-d <#,#,#>] [-M <trainP file>] [-h]\n",selfname);
	fprintf(stderr,"\tRequired:\n");
	fprintf(stderr,"\t-m <modelfile>\tmodel file, a line for mu and a line for pi, three\n");
	fprintf(stderr,"\t\t\tspace-separated values each, like:\n");
	fprintf(stderr,"\t\t\tmu 0.xxxxxxxx 0.xxxxxxxx 0.xxxxxxxx\n");
	fprintf(stderr,"\t\t\tpi 0.xxxxxxxx 0.xxxxxxxx 0.xxxxxxxx\n");
	fprintf(stderr,"\tOptional:\n");
	fprintf(stderr,"\t-i <infile>\tquality pileup, from pileup2matlab.pl script (def: STDIN)\n");
	fprintf(stderr,"\t-o <outfile>\twhere to put the results (def: STDOUT)\n");
	fprintf(stderr,"\t-T | -C | -F\tTrain, Classify or Filter (def: Classify)\n");
	fprintf(stderr,"\t-p q|m|s\tInput pileup format (def: s)\n\t\t\tq = quality\n\t\t\tm = MAQ\n\t\t\ts = SAMtools\n");
	fprintf(stderr,"\t-t mb|m|b|M|Mb|MB|SNVMix1\n\t\t\tFilter type (def: mb)\n");
	fprintf(stderr,"\t\t\tmb\tLowest between map and base quality\n");
	fprintf(stderr,"\t\t\tm\tFilter on map, and use as surrogate for base quality\n");
	fprintf(stderr,"\t\t\tb\tFilter on base quality, take map quality as 1\n");
	fprintf(stderr,"\t\t\tM\tFilter on map quality but use only base quality\n");
	fprintf(stderr,"\t\t\tMb\tFilter on map quality and use both map and base qualities\n");
	fprintf(stderr,"\t\t\tMB\tFilter on map quality AND base quality\n");
	fprintf(stderr,"\t\t\tSNVMix1\tFilter on base quality and map quality, afterwards consider them perfect\n");
	fprintf(stderr,"\t-q #\t\tCutoff Phred value for Base Quality, anything <= this value is ignored (def: Q%d)\n",default_params.bQ);
	fprintf(stderr,"\t-Q #\t\tCutoff Phred value for Map Quality, anything <= this value is ignored (def: Q%d)\n",default_params.mQ);
	fprintf(stderr,"\n\tTraining parameters:\n");
	fprintf(stderr,"\t-a #,#,#\tProvide alpha training parameters\n");
	fprintf(stderr,"\t-b #,#,#\tProvide beta training parameters\n");
	fprintf(stderr,"\t-d #,#,#\tProvide delta training parameters\n");
	fprintf(stderr,"\t-M <file>\tProvide a file containing training parameters (will override -a, -b and -d)\n");
	fprintf(stderr,"\t\t\tValues are space-separated:\n");
	fprintf(stderr,"\t\t\talpha # # #\n");
	fprintf(stderr,"\t\t\tbeta # # #\n");
	fprintf(stderr,"\t\t\tdelta # # #\n");
	fprintf(stderr,"\n\t-h\t\tthis message\n");
	fprintf(stderr,"\nImplementation: Rodrigo Goya, 2009. Version %s\n",VERSION);
	exit(0);
}

// Based on classify pileup, but reads the entire file to memory for training purposes.
void snvmixGetTrainSet_pileup(param_struct *params) {
//MAQ:
//	1       554484  C       1752    @,,.,.,, @xxx @xxx xxxxx
//SAM:
//	1       554484  C       1752    ,,.,.,,	xxx xxxx

	FILE *fptrIN;

	char *line = NULL;
	size_t line_len = 0, prev_len = 0;
	int read_len = 0;
	int col_num = 0;
	long int line_num = 0;

	char *col, *qual;
	char *col_sep = "\t", *col_str, *col_ptr;
	char *qual_sep = ",", *qual_str, *qual_ptr;
	char *chr, *pos, *ref, nref, call;	

	char *bQ, *mQ;
	int *calls, *tmpQ;
	int depth = 0, qual_num,  maxp;
	int depth_allocated = 0;

	int trainset = 0;
	int trainset_allocated = 0;

	int nref_count = 0, ref_count = 0;
	int i, call_i;

	fptrIN = params->input;

	if( (params->trainP.bQ = malloc(START_QNUM * sizeof (unsigned char **))) == NULL) {
		perror("ERROR allocating initial space for train.bQ"); exit(1);
	}
	if( (params->trainP.mQ = malloc(START_QNUM * sizeof (unsigned char **))) == NULL) {
		perror("ERROR allocating initial space for train.mQ"); exit(1);
	}
	if( (params->trainP.calls = malloc(START_QNUM * sizeof (signed char **))) == NULL) {
		perror("ERROR allocating initial space for train.calls"); exit(1);
	}
	if( (params->trainP.depth = malloc(START_QNUM * sizeof (int *))) == NULL) {
		perror("ERROR allocating initial space for train.depth"); exit(1);
	}
	trainset_allocated = START_QNUM;

	if( (calls = malloc(START_QNUM * sizeof (int))) == NULL) {
		perror("ERROR allocating initial space for train.depth"); exit(1);
	}
	depth_allocated = START_QNUM;
#if defined __linux__ || defined _GNU_SOURCE
	while( (read_len = getline(&line,&line_len,fptrIN)) != -1 ) {
#endif
#ifdef __APPLE__
	while( (line = fgetln(fptrIN,&line_len)) != NULL ) {
		line[line_len-1] = '\0';
#endif
		depth = 0;
		nref = 0;
		col_ptr = NULL;
		for(col_num = 0, col_str = line;  ; col_num++, col_str = NULL) {
			col = strtok_r(col_str, "\t", &col_ptr);
			if(col == NULL) {
				break;
			}
			if(col_num == 0) {
				chr = col;
			} else if(col_num == 1) {
				pos = col;
			} else if(col_num == 2) {
				ref = col;
			} else if(col_num == 3) {
				depth = atoi(col);
				if(depth > depth_allocated) {
					if( (calls = realloc(calls, sizeof (int) * (depth + START_QNUM))) == NULL) {
						perror("ERROR allocating space for calls"); exit(1);
					}
					depth_allocated = depth + START_QNUM;
				}
			} else if(col_num == 4) {
				if(params->input_type == M_PILEUP)
					col = col+sizeof(char);
					snvmixGetCallString(col, calls, depth, &nref);
			} else if(col_num == 5) {
				bQ = col;
				// If it's a MAQ pileup, we need to skip the @ sign
				if(params->input_type == M_PILEUP)
					bQ = bQ + sizeof(char);
				for(call_i = 0; call_i < depth; call_i++)
					bQ[call_i] = bQ[call_i]-33;
			} else if(col_num == 6) {
				mQ = col;
				// If it's a MAQ pileup, we need to skip the @ sign
				if(params->input_type == M_PILEUP)
					mQ = mQ + sizeof(char);
				for(call_i = 0; call_i < depth; call_i++)
					mQ[call_i] = mQ[call_i]-33;
			}
		}
		if(col_num != 7) {
			fprintf(stderr, "\nERROR: Pileup line %d has wrong number of columns (has %d expected 7, was flag '-s' used in samtools to create pileup?)\n", line_num, col_num);
			exit(1);
		}
		if(line_num >= trainset_allocated) {
			if( ( params->trainP.bQ = realloc(params->trainP.bQ, (line_num + START_QNUM) * sizeof (unsigned char *)) ) == NULL ) {
				perror("ERROR reallocating space for trainP.bQ"); exit(1);
			}
			if( ( params->trainP.mQ = realloc(params->trainP.mQ, (line_num + START_QNUM) * sizeof (unsigned char *)) ) == NULL ) {
				perror("ERROR reallocating space for trainP.mQ"); exit(1);
			}
			if( ( params->trainP.calls = realloc(params->trainP.calls, (line_num + START_QNUM) * sizeof (signed char *)) ) == NULL) {
				perror("ERROR reallocating space for trainP.calls"); exit(1);
			}
			if( ( params->trainP.depth = realloc(params->trainP.depth, (line_num + START_QNUM) * sizeof (int)) ) == NULL) {
				perror("ERROR reallocating space for trainP.depth"); exit(1);
			}
			trainset_allocated = line_num + START_QNUM;
		}

		nref = snvmixFilterCalls(calls,depth,bQ,mQ,params);
		nref_count = 0;
		ref_count = 0;
		for(qual_num = 0; qual_num < depth; qual_num++) {
			if(calls[qual_num] == -1)
				continue;
			if(calls[qual_num] == 1)
				ref_count++;
			if(calls[qual_num] == nref)
				nref_count++;
		}
		params->trainP.depth[line_num] = ref_count + nref_count;

		if( (params->trainP.bQ[line_num] = malloc(sizeof(unsigned char) * params->trainP.depth[line_num])) == NULL ) {
			perror("ERROR allocating space for trainP.bQ"); exit(1);
		}
		if( (params->trainP.mQ[line_num] = malloc(sizeof(unsigned char) * params->trainP.depth[line_num])) == NULL ) {
			perror("ERROR allocating space for trainP.mQ"); exit(1);
		}
		if( (params->trainP.calls[line_num] = malloc(sizeof(signed char) * params->trainP.depth[line_num])) == NULL ) {
			perror("ERROR allocating space for trainP.calls"); exit(1);
		}

		call_i = 0;
		for(qual_num = 0; qual_num < depth; qual_num++) {
			if(calls[qual_num] == -1)
				continue;
	
			params->trainP.calls[line_num][call_i] = (signed char) calls[qual_num];
			params->trainP.bQ[line_num][call_i] = (unsigned char) bQ[qual_num];
			params->trainP.mQ[line_num][call_i] = (unsigned char) mQ[qual_num];
			call_i++;
		}
		if( params->trainP.depth[line_num] != call_i ) {
			fprintf(stderr, "ERROR: mismatch between trainP.depth and call_i\n"); exit(1);
		}
		line_num++;
	}
	params->trainP.len = line_num;
	params->trainP.len = line_num;
	fclose(fptrIN);
	free(line);
}

// Use EM algorithm on a memory-located data set to train the parameters
void snvmixTrain_pileup(param_struct *params) {
	int line_num = 0, call_i = 0, k = 0;
	int iter = 0;
	FILE *fptrOUT;

	long double phred[PHRED_MAX + 1];
	long double bsum = 0, msum = 0, b[3], mu[3], notmu[3], pi[3], log_pi[3], z = 0;
	long double **pG, **xr;
	long double xrhat[3], sum_pG[3];
	long double Y, not_Y, M;
	int N_sum;
	int strength;

	long double ll, prev_ll = 0;
	//initPhred(phred, PHRED_MAX+1);
	initPhred(phred, PHRED_MAX);


	snvmixGetTrainSet_pileup(params);

	if(params->debug) {
	for(line_num = 0; line_num < params->trainP.len; line_num++) {
		fprintf(stderr, "line_num: %d", line_num);
		for(call_i = 0; call_i < params->trainP.depth[line_num]; call_i++) {
			fprintf(stderr, "\t(%d,bQ%d,mQ%d)", params->trainP.calls[line_num][call_i], params->trainP.bQ[line_num][call_i], params->trainP.mQ[line_num][call_i]);
		}
		fprintf(stderr, "\n\n");
	}
	}
	// Initialize mu and pi

	fprintf(stderr, "Initializing mu\n");
	for(k = 0; k < 3; k++) {
		params->mu[k] = (long double) params->trainP.alpha[k] / (params->trainP.alpha[k] + params->trainP.beta[k]);
		fprintf(stderr,"\talpha[%d] = %Lf,\tbeta[%d] = %Lf,\tmu[%d] = %Lf\n", k, params->trainP.alpha[k], k, params->trainP.beta[k], k, params->mu[k]);
	}

	fprintf(stderr, "Initializing pi\n");
	z = (long double) params->trainP.delta[0] + params->trainP.delta[1] + params->trainP.delta[2];
	if(!z) { z = 1; }
	for(k = 0; k < 3; k++) {
		params->pi[k] = (long double) params->trainP.delta[k]  / z;
		fprintf(stderr,"\tdelta[%d] = %Lf,\tpi[%d] = %Lf\n", k, params->trainP.delta[k], k, params->pi[k]);
	}

	strength = params->trainP.len;

	// Initialize matrices
	if( (pG = malloc(sizeof(long double *) * params->trainP.len)) == NULL) {
		perror("ERROR allocating initial space for pG"); exit(1);
	}
	if( (xr = malloc(sizeof(long double *) * params->trainP.len)) == NULL) {
		perror("ERROR allocating initial space for xr"); exit(1);
	}
	for(line_num = 0; line_num < params->trainP.len; line_num++) {
		// K = 3 cells for each line_num
		if( (pG[line_num] = malloc(sizeof(long double) * 3)) == NULL) {
			perror("ERROR allocating space for pG"); exit(1);
		}
		if( (xr[line_num] = malloc(sizeof(long double) * 3)) == NULL) {
			perror("ERROR allocating space for xr"); exit(1);
		}
	}

	for(iter = 0; iter < params->trainP.max_iter; iter++) {
		// Reset values for this iteration
		for(k = 0; k < 3; k++) {
			mu[k] = params->mu[k];
			notmu[k] = 1 - mu[k];
			pi[k] = params->pi[k];
			log_pi[k] = logl(pi[k]);

			sum_pG[k] = 0;
			xrhat[k] = 0;

		}
		N_sum = 0;
		ll = 0;

		// E step
		for(line_num = 0; line_num < params->trainP.len; line_num++) {
			if(params->trainP.depth == 0)
				continue;
			for(k = 0; k < 3; k++) {
				xr[line_num][k] = 0;
				b[k] = 0;
			}

			for(call_i = 0; call_i < params->trainP.depth[line_num]; call_i++) {
				if(params->trainP.calls[line_num][call_i] == 1) {
					Y = phred[params->trainP.bQ[line_num][call_i]];
					not_Y = 1-phred[params->trainP.bQ[line_num][call_i]];
				} else {
					Y = (1-phred[params->trainP.bQ[line_num][call_i]])/3;
					not_Y = phred[params->trainP.bQ[line_num][call_i]] + 2*(1-phred[params->trainP.bQ[line_num][call_i]])/3;
				}
				M =  phred[params->trainP.mQ[line_num][call_i]];
				for(k = 0; k < 3; k++) {
					b[k] = b[k] + logl( (1 - M) * 0.5 + M * (Y * mu[k] + not_Y * notmu[k]) );
					xr[line_num][k] = xr[line_num][k] + Y;
				}
			}
			z = 0;
			for(k = 0; k < 3; k++) {
				b[k] = expl(b[k] + log_pi[k]);
				z = z + b[k];
			}
			if(!z) { z=1; }
			for(k = 0; k < 3; k++) {
				pG[line_num][k] = b[k] / z;
			}

			ll = ll + logl(z);

			N_sum = N_sum + params->trainP.depth[line_num];

			for(k = 0; k < 3; k++) {
				sum_pG[k] = sum_pG[k] + pG[line_num][k];
				xrhat[k] = xrhat[k] + xr[line_num][k] * pG[line_num][k];
			}
		}

		// Check log likelihood
		if(iter == 0)
			prev_ll = ll;
		else if(ll <= prev_ll)
			break;
		prev_ll = ll;

		// M step
		z = 0;
		for(k = 0; k < 3; k++) {
			z = z + sum_pG[k] + params->trainP.delta[k];
		}
		if(!z) { z=1; }
		for(k = 0; k < 3; k++) {
			// Update pi
			params->pi[k] = (sum_pG[k] + params->trainP.delta[k]) / z;
			// Update mu
			params->mu[k] = (xrhat[k] + params->trainP.alpha[k]*strength-1) / (N_sum + params->trainP.alpha[k]*strength + params->trainP.beta[k]*strength-2);
		}
	}

	if(params->modelfile == NULL) {
		fptrOUT = stdout;
	} else {
		fptrOUT = fopen(params->modelfile, "w");
		if(fptrOUT == NULL) {
			perror("ERROR: could not open modelfile for writing, using STDOUT");
			fptrOUT = stdout;
		}
	}
	fprintf(fptrOUT,"mu");
	for(k = 0; k < 3; k++) {
		fprintf(fptrOUT, " %0.30Lf", params->mu[k]);
	}
	fprintf(fptrOUT,"\n");
	fprintf(fptrOUT,"pi");
	for(k = 0; k < 3; k++) {
		fprintf(fptrOUT, " %0.30Lf", params->pi[k]);
	}
	fprintf(fptrOUT,"\n");
}

