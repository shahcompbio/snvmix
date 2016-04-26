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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>

#define TYPE_mb 0
#define TYPE_m 1
#define TYPE_b 2
#define TYPE_M 3
#define TYPE_Mb 4
#define TYPE_MB 5
#define TYPE_SNVMix1 6

#define Q_PILEUP 0
#define M_PILEUP 1
#define S_PILEUP 2

#define PHRED_MAX 200

#define N 0
#define A 1
#define G 2
#define C 3
#define T 4

char base_code[] = {'N','A','G','C','T'};

struct params_train {
	long double alpha[3];
	long double beta[3];
	long double delta[3];
	char *param_file;
	int max_iter;
	unsigned char **bQ;
	unsigned char **mQ;
	signed char **calls;
	char *ref;
	int *pos;
	int *depth;
	int len;
};

typedef struct {
	FILE *input;
	FILE *output;
	char *inputfile;
	char *outputfile;
	char *modelfile;
	int filter_type;
	int train;
	int classify;
	int filter;
	int full;
	int input_type; // 0 = processed, 1 = maq pileup, 2 = sam pileup
	long double mu[3];
	long double pi[3];
	int max_iter;
	int bQ;
	int mQ;
	int debug;
	//struct {
	//	int alpha[3];
	//	int beta[3];
	//} train;
	struct params_train trainP;
} param_struct;

void updatePhred(long double *phredTable);
void initPhred(long double *phredTable, int elem);

void readTrainingParams(param_struct *params);
void readModel(param_struct *params);

void initSNVMix(int argc , char **argv, param_struct *params);
void resetParams(param_struct *params);

void usage(char *selfname);

void snvmixClassify_qualities(param_struct *params);
void snvmixClassify_pileup(param_struct *params);
void snvmixGetCallString(char *col, int *calls, int depth, char *nref);

int snvmixFilterCalls(int *calls, int depth, char *bQ, char *mQ, param_struct *params);
int snvmixSkipCall(int *calls, int qual_num, param_struct *params, char *bQ, char *mQ);

void snvmixTrain_qualities(param_struct *params);
void snvmixGetTrainSet_pileup(param_struct *params);
void snvmixTrain_pileup(param_struct *params);

long double normalise(long double *values, int len);
