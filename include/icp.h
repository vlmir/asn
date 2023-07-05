/*
iterative closest point
input: 2 sets of coordinates:
X (Besl:X) - target, reference, model...; immutable
Y (Besl:P) - source, reading, points...; rotated, translated
*/
#ifndef ICP_H
#define ICP_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h> // DBL_EPSILON
#include <math.h> // fabs, sqrt
#include <mpi.h>
#include "util.h"

typedef struct {
int rotmall[1]; // inirs
int twist[2]; // jrotsX, jrotsY
int puntaoptima[2]; // cpis, jrotsQ
int optimeproximo[64]; // x2z, frequency of using target pnts as closest ones
} Ret;

double RotMall(double **nefosX, double **nefosY, int ctrXpnts, int ctrYpnts,
               int wgtme, Ret *auxp);
int Twist(double RR[24][3][3], double CVX[3][3], double CVY[3][3], Ret *auxp);
double PuntaOptima(double **nuageX, double **nuageY, int ctrXpnts, int ctrYpnts,
                   int wgtme, Ret *auxp);
double *OptimeProximo(int d, int m, int n, double **x, double **y, double **z, Ret *auxp);

void rotate(int n, double S[n][n], int i, int j, int k, int l, double s,
            double tau) ;
int jacobi(int n, double S[n][n], double *w, double V[n][n],
           const int MAXITERS) ;
#endif
