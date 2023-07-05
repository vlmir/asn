#ifndef MPUTIL_H
#define MPUTIL_H

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

void PrintMp(char message[]);
void PrintLogMp(char *msg, const char *func);
void ErrorMp(char message[]);
void PrintC2dMp(char **src, int dim1, int dim2);
void PrintI2dMp(int **src, int dim1, int dim2);
int BcastD3d(double ***src, int dim1, int dim2, int dim3);
int BcastD2d(double **src, int dim1, int dim2);
int BcastI2d(int **src, int dim1, int dim2);
int OnelineD3d(double ***src, int i, int j);
int OnelineD2d(double **src, int i);
int ScattervI2d(int **global, int **local,      // source array, subarray
                const int ggs1, const int ggs2, // dimensions of global array
                const int lgs1, const int lgs2, // dimensions of subarray
                int root                        // ID of root process
                );
#endif
