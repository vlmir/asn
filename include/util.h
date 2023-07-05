#ifndef UTIL_H
#define UTIL_H

#define _GNU_SOURCE // asprintf
#include <stdlib.h>
#include <stdio.h>
#include <math.h> // sqrt
#include <stdarg.h>
#include <stddef.h>
#define NR_END 1
#define FREE_ARG char *

void NullD1d(int m, double *vec);
void NullD2d(int m, int n, double **arr);
void NullDm(int m, int n, double arr[m][n]);
void Error(char error_text[]);
void mpMpn(double **c, double **a, double **b, int m, int n, int p);
void PrintDm(int m, int n, double arr[m][n]);
void PrintDv(int m, double arr[m]);
void PrintD2d(double **arr, int m, int n);
void PrintI2d(int **arr, int m, int n);
void MeanDm(int m, int n, double **mat, double *vec);
void CopyDm(int m, int n, double src[m][n], double trg[m][n]);
void CopyD2d(int m, int n, double **src, double **trg);
int *Ivector(long nl, long nh);
int **Iarray2d(int m, int n);
int MallocC2d(char ***array, int dim1, int dim2);
char CallocC2d(char ***array, int dim1, int dim2);
int FreeC2d(char ***array);
int CallocD2d(double ***array, int dim1, int dim2);
int MallocD2d(double ***array, int dim1, int dim2);
int ReallocD2d(double ***array, int dim1, int dim2);
int FreeD2d(double ***array);
int MallocI2d(int ***array, int dim1, int dim2);
int CallocI2d(int ***array, int dim1, int dim2);
int FreeI2d(int ***array);
int CallocD3(double ****arr, long dim1, long dim2, long dim3);
int MallocD3(double ****arr, long dim1, long dim2, long dim3);
int FreeD3(double ****arr, long dim1, long dim2);
int MallocD3D(double ****arr, long dim1, long dim2, long dim3);
int FreeD3D(double ****arr, long dim1);
int MallocC3D(char ****arr, long dim1, long dim2, long dim3);
int FreeC3D(char ****arr, long dim1);
int MallocD3d(double ****arr, int dim1, int dim2, int dim3);
int CallocD3d(double ****arr, int dim1, int dim2, int dim3);
int FreeD3d(double ****arr);
int CmprDs(const void *a, const void *b);
int CmprIs(const void *a, const void *b);
int Count_inlines(const char *path);
int Read_table(char *path, char *fmtstr, int num, ...);
int **CentrifyD2d(double **object, int points, int dim);
int Squash(int points, int dim, double **obj, double arr[dim][dim]);
int Xwash(int n, int d, double **set1, double **set2, double xvar[d][d]);
int CountNLs(char *filepath);
int Chomp(char *str);
double *Dvector(long nl, long nh);
double **Dmatrix(long n1l, long n1h, long n2l, long n2h);
double **Darray2d(int m, int n); // vm: tested
double EDsquare(int dim, double *point1, double *point2);
double MSerror(int dim, int n, double **y, double **z, int w);
double Eudist(int dim, double *point1, double *point2);

void mulmat3(int nrow, int ncol, double a[nrow][ncol], int m, int n,
             double b[m][n], int nr, int nc, double c[nr][nc]) ;
double Determinant(int n, double a[n][n]) ;
void Transpose(int nrow, int ncol, double a[nrow][ncol], int m, int n,
               double b[m][n]) ;
int usage(const char *name, const char *help[], const int n) ;
#endif
