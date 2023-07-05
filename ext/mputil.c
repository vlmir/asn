/*
 * all the routines presume MPI running in the caller
 */
#include "mputil.h"

void PrintLogMp(char *msg, const char *func)
	//
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("%d:%s: %s %s: \n\t" 
			"%s\n",
      rank, func, __DATE__, __TIME__
			, msg);
  fflush(stdout);
}

void PrintMp(char message[]) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("%d: %s\n", rank, message);
  fflush(stdout);
}

void ErrorMp(char message[])
/* standard error handler */
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  fprintf(stderr, "MPI run-time error...\n");
  printf("rank%d: %s\n", rank, message);
  fprintf(stderr, "...now exiting ...\n");
  MPI_Abort(MPI_COMM_WORLD, 1);
}

void PrintC2dMp(char **src, int dim1, int dim2) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  for (int p = 0; p < size; p++) {
    if (rank == p) {
      printf("Local process on rank %d:\n", rank);
      fflush(stdout);
      for (int i = 0; i < dim1; i++) {
        putchar('|');
        for (int j = 0; j < dim2; j++) {
          putchar(src[i][j]);
        }
        printf("|\n");
        fflush(stdout);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void PrintI2dMp(int **src, int dim1, int dim2) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  for (int p = 0; p < size; p++) {
    if (rank == p) {
      printf("Local process on rank %d:\n", rank);
      fflush(stdout);
      for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
          printf("%2d ", src[i][j]);
          fflush(stdout);
        }
        printf("\n");
        fflush(stdout);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

int ScattervI2d(int **global, int **local,      // source array, subarray
                const int ggs1, const int ggs2, // dimensions of global array
                const int lgs1, const int lgs2, // dimensions of subarray
                int root                        // ID of root process
                ) {
  /* dimensions of grid of processes */
  const int pgs1 = ggs1 / lgs1;
  const int pgs2 = ggs2 / lgs2;

  /* sizes of arrays */
  int lgs = lgs1 * lgs2;
  int pgs = pgs1 * pgs2;

  /* rank of current process and no. of processes */
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (size != pgs) {
    fprintf(stderr, "Only works with np=%d for now\n", pgs);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* create a datatype to describe the subarrays of the global array */

  int globs[2] = {ggs1, ggs2}; /* size of global array */
  int locs[2] = {lgs1, lgs2};  /* size of sub-region */
  int starts[2] = {0, 0};      /* where this one starts */
  MPI_Datatype protype, subarr;
  MPI_Type_create_subarray(2, globs, locs, starts, MPI_ORDER_C, MPI_INT,
                           &protype);

  /* change the extent of the subarr */
  int ext; // the value is arbitrary, might be as well e.g.1
           // setting it to the length of subrows
  ext = lgs2;
  // NOTE for test files both pgs1 and lgs2 eq 3 - pure coincidence !!
  MPI_Type_create_resized(protype, 0, ext * sizeof(int), &subarr);
  MPI_Type_commit(&subarr);

  int *globalptr = NULL;
  if (rank == root)
    globalptr = &(global[0][0]);

  /* scatter the array to all processors */
  int counts[pgs]; // how many pieces of data everyone has
  // in units of blocks
  int displs[pgs]; // the starting point of everyone's data in the
                   // global array
                   // in block extents

  if (rank == root) {
    for (int i = 0; i < pgs; i++)
      counts[i] = 1;
    int disp = 0;
    for (int i = 0; i < pgs1; i++) {   // for each block in the 1st dimension
      for (int j = 0; j < pgs2; j++) { // for each block in the 2nd dimension
        displs[i * pgs2 + j] = disp;   // test: ind{0,1,2} val{0,5,10}
        disp += 1;
      }
      disp += (lgs1 - 1) * pgs2; // sic !!
    }
    printf("With displacements:\n");
    for (int i = 0; i < pgs1; i++) {
      for (int j = 0; j < pgs2; j++) {
        printf("%2d ", displs[i * pgs2 + j]);
      }
      printf("\n");
    }
  }

  MPI_Scatterv(               //
      globalptr,              //
      counts, displs, subarr, // proc i gets counts[i] types from displs[i]
      &(local[0][0]), lgs, MPI_INT, // receiving lgs MPI_INTs into local
      root, MPI_COMM_WORLD);

  /* now all processors print their local data: */
  PrintI2dMp(local, lgs1, lgs2);
  return 0;
}

int BcastD3d(double ***src, int dim1, int dim2, int dim3) {
  for (int i = 0; i < dim1; i++)
    for (int j = 0; j < dim2; j++)
      if (MPI_Bcast(src[i][j], dim3, MPI_DOUBLE, 0, MPI_COMM_WORLD) !=
          MPI_SUCCESS)
        return EXIT_FAILURE;
  return 0;
}

int BcastD2d(double **src, int dim1, int dim2) {
  for (int i = 0; i < dim1; i++)
    if (MPI_Bcast(src[i], dim2, MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
      return EXIT_FAILURE;
  return 0;
}

int BcastI2d(int **src, int dim1, int dim2) {
  for (int i = 0; i < dim1; i++)
    if (MPI_Bcast(src[i], dim2, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
      return EXIT_FAILURE;
  return 0;
}

int OnelineD3d(double ***src, int i, int j) {
  int comsize;
  if (MPI_Comm_size(MPI_COMM_WORLD, &comsize))
    comsize = 1;
  int myrank;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &myrank) != MPI_SUCCESS)
    myrank = 0;
  for (int r = 0; r < comsize; r++) {
    if (myrank == r) {
      printf("r:%d x:%lf y:%lf z:%lf \n", //
             myrank, src[i][j][0], src[i][j][1], src[i][j][2]);
    }
    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
      ErrorMp("MPI_Barrier failed"); // very rare case where barrier used
  }
  return 0;
}

int OnelineD2d(double **src, int i) {
  int comsize;
  if (MPI_Comm_size(MPI_COMM_WORLD, &comsize))
    comsize = 1;
  int myrank;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &myrank) != MPI_SUCCESS)
    myrank = 0;
  for (int r = 0; r < comsize; r++) {
    if (myrank == r) {
      printf("r:%d x:%lf y:%lf \n", //
             myrank, src[i][0], src[i][1]);
    }
    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
      ErrorMp("MPI_Barrier failed"); // very rare case where barrier used
  }
  return 0;
}

