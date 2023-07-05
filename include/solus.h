/*
 * Computing ALL pairwise distances between cloud points of two data sets
 * Currently the number of clouds in the sets is assimed identical
 */
#ifndef SOLUS_H
#define SOLUS_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <float.h> // DBL_EPSILON
#include <math.h>  // fabs, sqrt
#include <string.h>
#include <time.h>
//#include <unistd.h> // getopt?
#include <getopt.h> // or unistd.h?
#include "util.h"
#include "mputil.h"
#include "icp.h"

#define SYNOPSIS printf("synopsis: %s -f <file>\n", argv[0]) // TODO 
// all vars are presumed to be 'extern'
const int MASTER = 0; // default MPI master rank 
const int MAXPNTS = 32; // cutoff points/cloud used by strainer.pl 
/*
 * 'stout' controls STDOUT 
 * 0: STDOUT kept blank, only external messages
 * 1: plus summary info for solus 
 * 2: for debugging with VERY small test files 
*/
int stout = 1;

int Parse(FILE *fpr, double **data, int **meta, int mypnts);
int Spider(double **dataX, double **dataY, int **metaX, int **metaY, double **mydsts, int wgtme, int setctrs);
// currently the matrix is square 
//int WriteGrid(double **mydsts, char *binfn, int setctrs, int step);
#endif
