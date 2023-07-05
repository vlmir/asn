#include "util.h"

/// Multiplication of matrices for 2d array ///
void mulmat3(int nrow, int ncol, double a[nrow][ncol], int m, int n,
             double b[m][n], int nr, int nc, double c[nr][nc]) {
  int i, j, k;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < nc; j++)
      c[i][j] = 0.;
  }
  for (i = 0; i < nrow; i++) {
    for (j = 0; j < n; j++) {
      for (k = 0; k < m; k++) {
        c[i][j] += a[i][k] * b[k][j];
      }
    }
  }
}

////  Recursive definition of determinant using expansion by minors ////
double Determinant(int n, double a[n][n]) {
  int i, j, j1, j2;
  double det = 0;
  double m[n - 1][n - 1];
  if (n < 1) {
    // Error //
  } else if (n == 1) { // Shouldn't get used //
    det = a[0][0];
  } else if (n == 2) {
    det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  } else {
    det = 0;
    for (j1 = 0; j1 < n; j1++) {
      for (i = 1; i < n; i++) {
        j2 = 0;
        for (j = 0; j < n; j++) {
          if (j == j1)
            continue;
          m[i - 1][j2] = a[i][j];
          j2++;
        }
      }
      det += pow(-1.0, 1.0 + j1 + 1.0) * a[0][j1] * Determinant(n - 1, m);
    }
  }
  return (det);
}

////  Transpose of a square matrix, do it in situ ////
void Transpose(int nrow, int ncol, double a[nrow][ncol], int m, int n,
               double b[m][n]) {
  int i, j;
  for (i = 0; i < nrow; i++) {
    for (j = 0; j < ncol; j++) {
      b[j][i] = a[i][j];
    }
  }
}

// works as well with vec[] vars
void NullD1d(int m, double *vec) {
  for (int i = 0; i < m; i++)
    vec[i] = 0.;
}

// Note: doesn't accept arr[][]
void NullD2d(int m, int n, double **arr) {
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      arr[i][j] = 0.;
}
void NullDm(int m, int n, double arr[m][n]) {
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      arr[i][j] = 0.;
}

/*
 * essentially from NR
*/
void Error(char error_text[])
/* standard error handler */
{
  fprintf(stderr, "run-time error...\n");
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, "...now exiting to system...\n");
  exit(1);
}

// used only in playball.c
double *Dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v = (double *)calloc((size_t)(nh - nl + 1 + NR_END), (size_t)sizeof(double));
  if (!v)
    Error("allocation failure in Dvector()");
  return v - nl + NR_END;
}

// used only in pdbballs, 1 instance; to be replaced
int *Ivector(long nl, long nh)
/* allocate a int vector with subscript range v[nl..nh] */
{
  int *v;

  v = (int *)calloc((size_t)(nh - nl + 1 + NR_END), (size_t)sizeof(int));
  if (!v)
    Error("allocation failure in Ivector()");
  return v - nl + NR_END;
}

// used only in playball.c
double **Dmatrix(long n1l, long n1h, long n2l, long n2h)
/* allocate a double matrix with subscript range m[n1l..n1h][n2l..n2h] */
{
  long i, n11 = n1h - n1l + 1, n22 = n2h - n2l + 1;
  double **m;

  /* allocate pointers to rows */
  m = (double **)calloc((size_t)(n11 + NR_END), (size_t)sizeof(double *));
  if (!m)
    Error("allocation failure 1 in Dmatrix()");
  m += NR_END;
  m -= n1l;

  /* allocate rows and set pointers to them */
  m[n1l] =
      (double *)calloc((size_t)(n11 * n22 + NR_END), (size_t)sizeof(double));
  if (!m[n1l])
    Error("allocation failure 2 in Dmatrix()");
  m[n1l] += NR_END;
  m[n1l] -= n2l;

  for (i = n1l + 1; i <= n1h; i++)
    m[i] = m[i - 1] + n22;

  return m;
}

/*
 * the 3 subs below NOT USED
void Free_dvector(double *v, long nl)
// free a double vector allocated with Dvector() //
{
  free((FREE_ARG)(v + nl - NR_END));
}

void Free_ivector(int *v, long nl)
// free a double vector allocated with Dvector() //
{
  free((FREE_ARG)(v + nl - NR_END));
}

void Free_dmatrix(double **m, long nrl, long ncl)
// free a double matrix allocated by Dmatrix() //
{
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}
*/

/*
 * from
 * http://stackoverflow.com/questions/5585630/mpi-type-create-subarray-and-mpi-gather?answertab=votes#tab-top
 * of the same type as in NR
 * NOT USED
*/
int **Iarray2d(int m, int n) {
  int **arr = NULL;
  int i;

  arr = (int **)malloc(m * sizeof(int *));

  if (arr != NULL) {
    arr[0] = (int *)malloc(m * n * sizeof(int));
    if (arr[0] != NULL)
      for (i = 1; i < m; i++)
        arr[i] = arr[0] + i * n;
    else {
      free(arr);
      arr = NULL;
    }
  } else {
    arr = NULL;
  }
  return arr;
}

double **Darray2d(int m, int n) { // vm: tested
  double **arr = NULL;
  int i;

  arr = (double **)malloc(m * sizeof(double *));

  if (arr != NULL) {
    arr[0] = (double *)malloc(m * n * sizeof(double));
    if (arr[0] != NULL)
      for (i = 1; i < m; i++)
        arr[i] = arr[0] + i * n;
    else {
      free(arr);
      arr = NULL;
    }
  } else {
    arr = NULL;
  }
  return arr;
}

/*
 * from
 * http://stackoverflow.com/questions/9269399/sending-blocks-of-2d-array-in-c-using-mpi
 */
int MallocC2d(char ***array, int dim1, int dim2) {
  /* allocate the dim1*dim2 contiguous items */
  char *p = (char *)malloc(dim1 * dim2 * sizeof(char));
  if (!p)
    Error("MallocC2d:");

  /* allocate the row pointers into the memory */
  (*array) = (char **)malloc(dim1 * sizeof(char *));
  if (!(*array)) {
    free(p);
    Error("MallocC2d:");
  }

  /* set up the pointers into the contiguous memory */
  for (int i = 0; i < dim1; i++)
    (*array)[i] = &(p[i * dim2]);

  return 0;
}

char CallocC2d(char ***array, int dim1, int dim2) {
  int i;
  /* allocate the dim1*dim2 contiguous items */
  char *p = calloc(dim1 * dim2, sizeof(char));
  if (!p)
    Error("CallocC2d:");

  /* allocate the row pointers into the memory */
  (*array) = calloc(dim1, sizeof(char *));
  if (!(*array)) {
    free(p);
    Error("CallocC2d:");
  }

  /* set up the pointers into the contiguous memory */
  for (i = 0; i < dim1; i++)
    (*array)[i] = &(p[i * dim2]);

  return 0;
}

int FreeC2d(char ***array) {
  /* free the memory - the first element of the array is at the start */
  free(&((*array)[0][0]));
  /* free the pointers into the memory */
  free(*array);
  return 0;
}

int CallocD2d(double ***array, int dim1, int dim2) {
  // tested
  /* allocate the dim1*dim2 contiguous items */
  double *p = (double *)calloc(dim1 * dim2, sizeof(double));
  if (!p)
    Error("CallocD2d:1");

  /* allocate the row pointers into the memory */
  (*array) = calloc(dim1, sizeof(double *));
  if (!(*array)) {
    free(p);
    Error("CallocD2d:2");
  }

  /* set up the pointers into the contiguous memory */
  for (int i = 0; i < dim1; i++)
    (*array)[i] = &(p[i * dim2]);

  return 0;
}

int MallocD2d(double ***array, int dim1, int dim2) {
  // tested
  /* allocate the dim1*dim2 contiguous items */
  double *p = (double *)malloc(dim1 * dim2 * sizeof(double));
  if (!p)
    Error("MallocD2d:");

  /* allocate the row pointers into the memory */
  (*array) = malloc(dim1 * sizeof(double *));
  if (!(*array)) {
    free(p);
    Error("MallocD2d:");
  }

  /* set up the pointers into the contiguous memory */
  for (int i = 0; i < dim1; i++)
    (*array)[i] = &(p[i * dim2]);

  return 0;
}

int ReallocD2d(double ***array, int dim1, int dim2) {
  // seg fault
  /* reallocate the dim1*dim2 contiguous items */
  double *p = (double *)realloc((**array), dim1 * dim2 * sizeof(double));
  if (!p)
    Error("ReallocD2d:");

  /* allocate the row pointers into the memory */
  (*array) = realloc((*array), dim1 * sizeof(double *));
  if (!(**array)) {
    free(p);
    Error("ReallocD2d:");
  }

  /* set up the pointers into the contiguous memory */
  for (int i = 0; i < dim1; i++)
    (*array)[i] = &(p[i * dim2]);

  return 0;
}

int FreeD2d(double ***array) {
  /* free the memory - the first element of the array is at the start */
  free(&((*array)[0][0]));
  /* free the pointers into the memory */
  free(*array);
  return 0;
}

int MallocI2d(int ***array, int dim1, int dim2) {
  int i;
  /* allocate the dim1*dim2 contiguous items */
  int *p = malloc(dim1 * dim2 * sizeof(int));
  if (!p)
    Error("MallocI2d:");

  /* allocate the row pointers into the memory */
  (*array) = malloc(dim1 * sizeof(int *));
  if (!(*array)) {
    free(p);
    Error("MallocI2d:");
  }

  /* set up the pointers into the contiguous memory */
  for (i = 0; i < dim1; i++)
    (*array)[i] = &(p[i * dim2]);

  return 0;
}

int CallocI2d(int ***array, int dim1, int dim2) {
  int i;
  /* allocate the dim1*dim2 contiguous items */
  int *p = calloc(dim1 * dim2, sizeof(int));
  if (!p)
    Error("CallocI2d:");

  /* allocate the row pointers into the memory */
  (*array) = calloc(dim1, sizeof(int *));
  if (!(*array)) {
    free(p);
    Error("CallocI2d:");
  }

  /* set up the pointers into the contiguous memory */
  for (i = 0; i < dim1; i++)
    (*array)[i] = &(p[i * dim2]);

  return 0;
}

int FreeI2d(int ***array) {
  /* free the memory - the first element of the array is at the start */
  free(&((*array)[0][0]));
  /* free the pointers into the memory */
  free(*array);
  return 0;
}

/*----------------------- MINE ----------------------------------------------*/

int CallocD3(double ****arr, long dim1, long dim2, long dim3)
/*
 * allocate a 3-dimentional structure for doubles
 * subscript range dim1[0..dim1-1][0..dim2-1][0..dim3-1]
*/
{
  (*arr) = (double ***)calloc(dim1, sizeof(double **));
  if (!(*arr))
    Error("CallocD3: failure dimension 1");
  for (long i = 0; i < dim1; i++) {
    (*arr)[i] = (double **)calloc(dim2, sizeof(double *));
    if (!(*arr)[i])
      Error("CallocD3: failure dimension 2");
    for (long j = 0; j < dim2; j++) {
      (*arr)[i][j] = (double *)calloc(dim3, sizeof(double));
      if (!(*arr)[i][j])
        Error("CallocD3: failure dimension 3");
    }
  }
  return 0;
}

int MallocD3(double ****arr, long dim1, long dim2, long dim3)
/*
 * allocate a 3-dimentional structure for doubles
 * subscript range dim1[0..dim1-1][0..dim2-1][0..dim3-1]
*/
{
  (*arr) = (double ***)malloc(dim1 * sizeof(double **));
  if (!(*arr))
    Error("MallocD3: failure dimension 1");
  for (long i = 0; i < dim1; i++) {
    (*arr)[i] = (double **)malloc(dim2 * sizeof(double *));
    if (!(*arr)[i])
      Error("MallocD3: failure dimension 2");
    for (long j = 0; j < dim2; j++) {
      (*arr)[i][j] = (double *)malloc(dim3 * sizeof(double));
      if (!(*arr)[i][j])
        Error("MallocD3: failure dimension 3");
    }
  }
  return 0;
}

int FreeD3(double ****arr, long dim1, long dim2) {
  // NOTE: the order of statements if essential !
  for (long i = 0; i < dim1; i++)
    for (long j = 0; j < dim2; j++)
      free((*arr)[i][j]);
  for (long i = 0; i < dim1; i++)
    free((*arr)[i]);
  free((*arr));
  return 0;
}

int MallocD3D(double ****arr, long dim1, long dim2, long dim3)
/*
 * allocate dim1 blocks of contiguous dim2*dim3 items
 * subscript range dim1[0..dim1-1][0..dim2-1][0..dim3-1]
*/
{
  (*arr) = (double ***)malloc(dim1 * sizeof(double **));
  if (!(*arr))
    Error("MallocD3D: failure dimension 1");
  for (long i = 0; i < dim1; i++) {
    double **a;
    MallocD2d(&a, dim2, dim3);
    (*arr)[i] = a;
  }
  return 0;
}

int FreeD3D(double ****arr, long dim1) {
  // NOTE: the order of statements if essential !
  for (long i = 0; i < dim1; i++)
    free((*arr)[i]);
  free((*arr));
  return 0;
}

int MallocC3D(char ****arr, long dim1, long dim2, long dim3)
/*
 * allocate dim1 blocks of contiguous dim2*dim3 items
 * subscript range dim1[0..dim1-1][0..dim2-1][0..dim3-1]
*/
{
  (*arr) = (char ***)malloc(dim1 * sizeof(char **));
  if (!(*arr))
    Error("MallocC3D: failure dimension 1");
  for (long i = 0; i < dim1; i++) {
    char **a;
    MallocC2d(&a, dim2, dim3);
    (*arr)[i] = a;
  }
  return 0;
}

int FreeC3D(char ****arr, long dim1) {
  // NOTE: the order of statements if essential !
  for (long i = 0; i < dim1; i++)
    free((*arr)[i]);
  free((*arr));
  return 0;
}

int MallocD3d(double ****arr, int dim1, int dim2, int dim3) {
  /* allocate the dim1*dim2*dim3 contiguous items */
  // 2016-11-19: works but the output different from MallocD3D
  double *p = (double *)malloc(dim1 * dim2 * dim3 * sizeof(double));
  if (!p)
    Error("MallocD3d: failure allocating contiguous memory");
  (*arr) = (double ***)malloc(dim1 * sizeof(double **));
  if (!(*arr)) {
    free(p);
    Error("MallocD3d: failure dimension 1");
  }
  for (int i = 0; i < dim1; i++) {
    (*arr)[i] = (double **)malloc(dim2 * sizeof(double *));
    if (!(*arr)[i]) {
      free(p);
      Error("MallocD3d: failure dimension 2");
    }
    for (int j = 0; j < dim2; j++) {
      (*arr)[i][j] = &(p[(i * dim3 + j) * dim3]); // TODO: test pointer location
    }
  }
  return 0;
}

int CallocD3d(double ****arr, int dim1, int dim2, int dim3) {
  /* allocate the dim1*dim2*dim3 contiguous items */
  double *p = (double *)calloc(dim1 * dim2 * dim3, sizeof(double));
  if (!p)
    Error("CallocD3d: failure allocating contiguous memory");
  (*arr) = (double ***)calloc(dim1, sizeof(double **));
  if (!(*arr)) {
    free(p);
    Error("CallocD3d: failure dimension 1");
  }
  for (int i = 0; i < dim1; i++) {
    (*arr)[i] = (double **)calloc(dim2, sizeof(double *));
    if (!(*arr)[i]) {
      free(p);
      Error("CallocD3d: failure dimension 2");
    }
    for (int j = 0; j < dim2; j++) {
      (*arr)[i][j] = &(p[(i * dim3 + j) * dim3]); // TODO: test pointer location
    }
  }
  return 0;
}

int FreeD3d(double ****arr) {
  // tested 2016-10-25
  /* free the memory - the first element of the array is at the start */
  free(&((*arr)[0][0][0]));
  /* free the pointers into the memory */
  free(*arr);
  return 0;
}
// to be used by qsort
int CmprDs(const void *a, const void *b) {
  if (*(double *)a < *(double *)b)
    return -1;
  if (*(double *)a > *(double *)b)
    return 1;
  if (*(double *)a == *(double *)b)
    return 0;
  return 2;
}

// to be used by qsort
int CmprIs(const void *a, const void *b) {
  if (*(int *)a < *(int *)b)
    return -1;
  if (*(int *)a > *(int *)b)
    return 1;
  if (*(int *)a == *(int *)b)
    return 0;
  return 2;
}

int Count_inlines(const char *path) {
  FILE *fp;
  if ((fp = fopen(path, "r")) == NULL) {
    fprintf(stderr, "Failed to read file %s\n", path);
    return EXIT_FAILURE;
  }
  int i = 0;
  int ch = 0;
  while ((ch = getc(fp)) != EOF) {
    if (ch == '\n')
      i++;
  }
  fclose(fp);
  return i;
}

int Read_table(char *path, char *fmtstr, int num, ...) {
  // a draft
  va_list valist;
  /* initialize valist for num number of arguments */
  va_start(valist, num); // 'num' - the argument before the ellipsis
  /* access all the arguments assigned to valist */
  void *args[num];
  int i = 0;
  for (i = 0; i < num; i++) {
    args[i] = va_arg(valist, void *);
  }

  FILE *fp;
  if ((fp = fopen(path, "r")) == NULL) {
    fprintf(stderr, "Failed to read file %s\n", path);
    return EXIT_FAILURE;
  }
  printf("Reading file %s\n", path);
  while (!feof(fp)) {
    int j;
    for (j = 0; j < num; j++) {
      if ((fscanf(fp, fmtstr, args)) < num) {
        fprintf(stderr, "failed assignments in fscanf");
        i++;
      }
    }
  }
  fclose(fp);
  /* clean memory reserved for valist */
  va_end(valist);
  return 0;
}

/* Alternative version
int **CentrifyD2d(int points, int dim, double **subject, double **object, int
ego) {
        for (int i=0; i<points; i++) for (int j=0; j<dim; j++)
object[i][j]=subject[i][j];
  if (ego < 0) { // centrify round the mean
    for (int k = 0; k < dim; k++) {
      double sum = 0.;
      for (int i = 0; i < points; i++)
        sum += object[i][k];
      for (int i = 0; i < points; i++)
        object[i][k] = object[i][k] - sum / (double)points;
    }
  } else { // centrify round the point with the index ego 0..points
    double self[dim];
    for (int k = 0; k < dim; k++)
      self[k] = object[ego][k];
    for (int i = 0; i < points; i++)
      for (int k = 0; k < dim; k++)
        object[i][k] -= self[k];
  }
  return 0;
}
*/

// centrify round the mean
int **CentrifyD2d(double **object, int points, int dim) {
    for (int k = 0; k < dim; k++) {
      double sum = 0.;
      for (int i = 0; i < points; i++)
        sum += object[i][k];
      for (int i = 0; i < points; i++)
        object[i][k] = object[i][k] - sum / (double)points;
  }
  return 0;
}

/// Covariance matrix <- double pointer ///
int Squash(int points, int dim, double **obj, double arr[dim][dim]) {
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      arr[i][j] = 0.;
      for (int row = 0; row < points; row++) {
        arr[i][j] += obj[row][i] * obj[row][j];
      }
    }
  }
  return 0;
}

int Xwash(int n, int d, double **set1, double **set2, double xvar[d][d]) {
  for (int k = 0; k < d; k++)
    for (int l = 0; l < d; l++)
      xvar[k][l] = 0.;
  for (int k = 0; k < d; k++)
    for (int l = 0; l < d; l++)
      for (int j = 0; j < n; j++)
        xvar[k][l] += set1[j][k] * set2[j][l];
  for (int k = 0; k < d; k++) {
    for (int l = 0; l < d; l++)
      xvar[k][l] = xvar[k][l] / (double)n;
  }
  return 0;
}

// count lines in file(s)
int CountNLs(char *filepath) {
  /* if filepath is a glob (e.g. *.xyz) returns total count */
  FILE *fp;
  char *command;
  if (asprintf(&command, "wc -l %s", filepath) == -1)
    Error("CountNLs:");
  /* Open the command for reading. */
  fp = popen(command, "r");
  if (fp == NULL) {
    printf("Failed to run command\n");
    exit(1);
  }
  /* Read the output a line at a time */
  int count = 0; // number of lines
  char *string;  // file path or 'total'
  if ((string = (char *)calloc(128, sizeof(char))) == NULL)
    Error("CountNLs: calloc failed");
  while (!feof(fp))
    if (fscanf(fp, "%d %s\n", &count, string))
      continue;
  pclose(fp);
  return count;
}

// square of Euclidean distance for 2 dim-dimensional points
double EDsquare(int dim, double *point1, double *point2) {
  double devi = .0;
  for (int k = 0; k < dim; k++)
    devi += (point1[k] - point2[k]) * (point1[k] - point2[k]);
  return devi;
}

/*
 * y, z: arrays of n dim-dimensional points
 * w: weighting factor, 0: no weighting, n: all pairs of points weighted
 */
/// computes Mean Square Error between source and target (RMSD sqared)
double MSerror(int dim, int n, double **y, double **z, int w) {
  double mse;         // mean sum of squared distances, returned
  double ssd = 0.;    // sum of squared distances
  double *y2z = NULL; // squared dstances between y and z points
  if ((y2z = (double *)calloc(n, sizeof(double))) == NULL)
    Error("calloc failed");
  for (int i = 0; i < n; i++)
    y2z[i] = EDsquare(dim, y[i], z[i]);
  double ssw = 0.; // sum of squared weights// TODO sort out
  if (!w) {
    for (int i = 0; i < n; i++)
      ssd += y2z[i];
    ssw = n; // ssw==n if all wgt==1
  } else {
    double wgt = 0.;
    qsort(y2z, n, sizeof(double), CmprDs);
    // split n in perikentrikos + apokentrikos:
    int perikentrikos = n - w;
    int apokentrikos = w;
    // smaller mindsts unweighted
    for (int i = 0; i < perikentrikos; i++) {
      ssd += y2z[i];
      ssw += 1; // weights set to '1'
    }
    // larger mindists weighted
    // weight changing from 1-1/apokentrikos to 0
    /* ori:
    for (int i = perikentrikos; i < n; i++) { // ori
    wgt = ((n - i - 1) * (n - i - 1)) / (apokentrikos * apokentrikos);
    */
    int myshift = 0; // 0: identical to ori:
    for (int i = 0; i < apokentrikos; i++) {
      wgt = 1 - (i + 1) / (apokentrikos + myshift);
      /* ori:
      wgt = (i+1)/apokentrikos -1 // effectively eliminates the
      last point
      w==n? sqrt(wgt) ~ (1-1/n)..0
      */
      ssd += y2z[i + perikentrikos] * wgt * wgt;
      ssw += wgt * wgt;
    }
  }
  mse = ssd / (double)ssw;
  return mse;
}
///////////////////////////////////////////////////////////////////////////
double Eudist(int dim, double *point1, double *point2) {
  double devi = .0;
  for (int k = 0; k < dim; k++)
    devi += (point1[k] - point2[k]) * (point1[k] - point2[k]);
  return sqrt(devi);
}

int Chomp(char *str) {
  /* replaces the first \n in the input string with \0 */
  int i;
  int flag = 0;
  int length = sizeof(str) + 1;
  for (i = 0; i < length; i++) {
    if (str[i] == '\n') {
      str[i] = '\0';
      flag++;
      break;
      return 0;
    }
  }
  if (flag == 0)
    fprintf(stderr, "No EOL in string '%s'\n", str);
  return 1;
}

void mpMpn(double **c, double **a, double **b, int m, int n, int p) {
  // to be tested TODO
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      c[i][j] = 0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++)
        c[i][j] += a[i][k] * b[k][j];
}

void PrintDm(int m, int n, double arr[m][n]) {
  printf("matrix of doubles size:%d*%d \n", m, n);
  fflush(stdout);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      printf("%lf ", arr[i][j]);
      fflush(stdout);
    }
    printf("\n");
    fflush(stdout);
  }
}

void PrintDv(int m, double arr[m]) {
  printf("vector of doubles size:%d \n", m);
  fflush(stdout);
  for (int i = 0; i < m; i++) {
    printf("%lf ", arr[i]);
    fflush(stdout);
  }
  printf("\n");
  fflush(stdout);
}

void PrintSv(int m, char * arr[m]) {
// TODO to be tested
  printf("vector of strings size:%d \n", m);
  fflush(stdout);
  for (int i = 0; i < m; i++) {
    printf("%s ", arr[i]);
    fflush(stdout);
  }
  printf("\n");
  fflush(stdout);
}

// Note: doesn't accept arr[][]
void PrintD2d(double **arr, int m, int n) {
  printf("matrix of doubles size:%d*%d \n", m, n);
  fflush(stdout);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      printf("(%d,%d) %1.3lf ", i, j, arr[i][j]);
      fflush(stdout);
    }
    printf("\n");
    fflush(stdout);
  }
}

// Note: doesn't accept arr[][]
void PrintI2d(int **arr, int m, int n) {
  printf("matrix of ints size:%d*%d \n", m, n);
  fflush(stdout);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      printf("%d ", arr[i][j]);
      fflush(stdout);
    }
    printf("\n");
    fflush(stdout);
  }
}

void MeanDm(int m, int n, double **mat, double *vec) {
  for (int j = 0; j < n; j++) {
    vec[j] = 0.;
    for (int i = 0; i < m; i++)
      vec[j] += mat[i][j];
    vec[j] /= (double)m;
  }
}

void CopyDm(int m, int n, double src[m][n], double trg[m][n]) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      trg[i][j] = src[i][j];
  }
}

void CopyD2d(int m, int n, double **src, double **trg) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      trg[i][j] = src[i][j];
  }
}

//int usage(const char *name, const char **help, const int n)
int usage(const char *name, const char *help[], const int n) // the same but clearer
{
	printf("Usage: ");
	printf("%s [-afhv -d <string>] -r <int> -c <int> -i <sting>\n", name);

	for (int i=0; i<n; i++)
	{
		printf("%s\n", help[i]);
	};
	return (0);
}
