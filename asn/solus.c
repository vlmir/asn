#include "solus.h"

int Parse(FILE *fpr, double **data, int **meta, int xclu)
// 
{

  int i = -1; // ctr index
  int minpnts = 2;
  int j = minpnts - 1; // pnt index within a center
  int l = 0;           // line index in the OUTPUT
  int xd = -1;
  /// vars for reading:
  int w;             // set-wide center id
  double x, y, z, d; // coordinates,  distance from the center
  int t;             // center atom type index
  int v;             // entry-wide atom index (original)
  int u;             // set-wide entry index, currently not used
  char *fmt = "%d %lf %lf %lf %lf %d %d %d\n";
  while (!feof(fpr))
  //	for (int l = 0; feof(fpr); l++) // didn't work = why??
  {
    if ((fscanf(fpr, fmt, &w, &x, &y, &z, &d, &t, &v, &u)) < 8)
      Error("failed fscanf ");
    xd = ceil(d); // == 0 for the center atom ONLY
    if (stout > 1)
      printf(">> v:%d i:%d j:%d l:%d xd:%d \n", v, i, j, l, xd);
    if (xd == 0) // || (w>i); new ctr
    {
      if (j > minpnts - 2) //
      {
        i++;
        meta[i][0] = t;
        meta[i][1] = v;
        meta[i][2] = l; // set wide index of the center atom
      }
      j = 0; // resetting pnt counter
    }
    if (j < xclu) //
    {
      data[l][0] = x;
      data[l][1] = y;
      data[l][2] = z;
      data[l][3] = d;
      meta[i][3] = j + 1; // cloud points; used in Spider!!
      meta[i][4] = xd;    // for computing distribution
      l++;
    }
    if (stout > 1)
      printf("<< v:%d i:%d j:%d l:%d xd:%d \n", v, i, j, l, xd);
    j++;
  }
  int setctrs = i + 1;
  int setpnts = l; // retained
  if (stout) PrintLogMp("data loaded:", __func__);
  if (stout) printf("\tsetctrs:%d setpnts:%d \n", setctrs, setpnts);
  return setpnts;
}

///////////////////////////////////////////////////////////////////////////////

int Spider(double **dataX, double **dataY, int **metaX, int **metaY, //
             double **mydsts, int wgtme, int setctrs) 
{ // Generic network computatiom for 2 sets of point clouds 
  const int D = 3; // dimentionality, xyz in this case; TODO: as arg
  const int MAXITERS = 512;
  int rank;
  int pool;
  MPI_Comm_size(MPI_COMM_WORLD, &pool);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int step = setctrs % pool ? setctrs / pool + 1 : setctrs / pool;
  int mystart = rank * step;
  int myend = mystart + step; // Note: myend is non-inclusive
  if (myend > setctrs)
    myend = setctrs;
  int inirsmin = 24;
  int cpismin = MAXITERS;
  int inirsmax = 0;
  int cpismax = 0;
  int inirsall = 0;
  int cpisall = 0;

  // computed in RotMall:
  int inirs; // number of R matrices used effectively to reach dstmin

  // computed in PuntaOptima:
  int cpis;   // number of closest point iterations
  int jrotsQ; // rotations in jacobi, cumulative

  // computed in Twist:
  int jrotsX; // rotations in jacobi, cumulative
  int jrotsY; // rotations in jacobi, cumulative

  Ret aux;
  aux.rotmall[0] = 0;     // inirs
  aux.twist[0] = 0;       // jrotsX, jrotsY
  aux.twist[1] = 0;       // jrotsX, jrotsY
  aux.puntaoptima[0] = 0; // cpis, jrotsQ
  aux.puntaoptima[1] = 0; // cpis, jrotsQ

  double **nefosX;
  for (int u = mystart; u < myend; u++) {
    int ctrXpnts = metaX[u][3]; // map1[u][1]; TODO: use metaX instead
    int ctrXind = metaX[u][2];  // int poffset = u * ctrXpnts;
    CallocD2d(&nefosX, ctrXpnts, D);
    for (int i = 0; i < ctrXpnts; i++)
      for (int k = 0; k < D; k++)
        nefosX[i][k] = dataX[ctrXind + i][k];
    // pool==2: nefosX properly split

    CentrifyD2d(nefosX, ctrXpnts, D);
    double **nefosY;
    for (int v = 0; v < setctrs; v++) {
      int ctrYpnts = metaY[v][3];
      int ctrYind = metaY[v][2];
      CallocD2d(&nefosY, ctrYpnts, D);
      for (int i = 0; i < ctrYpnts; i++)
        for (int k = 0; k < D; k++)
          nefosY[i][k] = dataY[ctrYind + i][k];
      /*
      if (stout) printf("u:%d v:%d ctrXind:%d ctrYind:%d \n", u, v, ctrXind, ctrYind );
      */
      CentrifyD2d(nefosY, ctrYpnts, D);
      // pool==2: nefosY properly split
      /////////////////////////////////////////////////////////////////////////
      //
      // RotMall calls PuntaOptima calls Twist
      double dstmin = RotMall(nefosX, nefosY, ctrXpnts, ctrYpnts, wgtme, &aux);
      //
      /////////////////////////////////////////////////////////////////////////
      /*
      if (stout) printf("dstmin:%f \n", dstmin);
      */
      mydsts[u - mystart][v] = dstmin;
      inirs = aux.rotmall[0];
      cpis = aux.puntaoptima[0];
      if (inirs < inirsmin)
        inirsmin = inirs;
      if (inirs > inirsmax)
        inirsmax = inirs;
      if (cpis < cpismin)
        cpismin = cpis;
      if (cpis > cpismax)
        cpismax = cpis;
      inirsall += inirs;
      cpisall += cpis;
    } // v<setctrs
    FreeD2d(&nefosY);
  } // u<myctrs
  jrotsQ = aux.puntaoptima[1];
  jrotsX = aux.twist[0];
  jrotsY = aux.twist[1];

  /// inirs: # R mtx used; cpis: # registrations in PuntaOptima
  if (stout) printf(
      "%d:%s: %s %s: \n\t"
      "inirsmin:%d inirsmax:%d inirsall:%d cpismin:%d cpismax:%d cpisall:%d \n",
      rank, __func__, __DATE__, __TIME__ //
      ,
      inirsmin, inirsmax, inirsall, cpismin, cpismax, cpisall);
  if (stout) printf("%d:%s: %s %s: \n\t"
         "jrotsX:%d jrotsY:%d jrotsQ:%d \n",
         rank, __func__, __DATE__, __TIME__ //
         ,
         jrotsX, jrotsY, jrotsQ);
  fflush(stdout);
  FreeD2d(&nefosX);
  FreeD2d(&dataX);
  // FreeD2d(&dataY);
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
// ICP 4 ASN - a single input file of multiple clouds 
{
  MPI_Init(&argc, &argv);
  int rank;
  int pool;
  MPI_Comm_size(MPI_COMM_WORLD, &pool);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int pariah = pool - 1;
  if (rank == pariah)
    if (stout) PrintLogMp("Let's start...", __func__);

  /// defaults for command line args
  char *infn = NULL;
  char *tabfn = NULL;
  int xclu = MAXPNTS;  // # of points per cloud to use; 2..MAXPNTS
  int wgtme = 0; // # of distances to weight, 0..MAXPNTS;
  int ictrs = 0; // # of clouds in the input
  int ierr = 0;

  extern char *optarg;
  int opt;
  if (rank == MASTER) {
    // the order of options is irrelevant
    char *opts = "i:o:c:x:w:h";
    while ((opt = getopt(argc, argv, opts)) != EOF)
      switch (opt) {
      case 'i':
        infn = optarg;
        break;
      case 'o':
        tabfn = optarg;
        break;
      case 'c':
        ictrs = atoi(optarg);
        break;
      case 'x':
        xclu = atoi(optarg); // clang: use 'strtol' if str > 1 char
        break;
      case 'w':
        wgtme = atoi(optarg);
        break;
      case 'h':
        ierr = 1;
        break;
      case '?':
        ierr = 1;
        break;
      }
    if (infn == NULL)
      ierr = 1;
    if (tabfn == NULL)
      ierr = 1;
    if (ictrs < 2)
      ierr = 1;
    if (xclu < 2)
      ierr = 1;
    if (ierr)
      SYNOPSIS;
    if (wgtme < 0)
      wgtme = xclu + wgtme + 1;
  }

  MPI_Bcast(&ierr, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&xclu, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&wgtme, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&ictrs, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  /* allocating memory for input data */
  int setctrs = ictrs;           // only one data set
  int setpnts = ictrs * MAXPNTS; // conservative estimate for Parse

  double **data;
  CallocD2d(&data, setpnts, 4);
  int **meta;
  CallocI2d(&meta, setctrs, 5);

  /* loading data */
  if (rank == MASTER) {
    FILE *fpr;
    if ((fpr = fopen(infn, "r")) == NULL)
      fprintf(stderr, "%s:fpr: failed to open file: %s \n", __func__, infn);

    ///////////////////////////////////////////////////////////////////////////
    setpnts = Parse(fpr, data, meta, xclu); // real value of setpnts
    ///////////////////////////////////////////////////////////////////////////

    if (!setpnts)
      Error("No data in the input");
    fclose(fpr);
    /// No barrier ! // why?
  } // Read data on rank 0

  /*
  * Sharing data from rank 0 (executed by EACH rank)
  * each process needs FULL data set;
  * empty lines at the bottom of data excluded
  */
  if (MPI_Bcast(&setctrs, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Bcast failed");
  if (MPI_Bcast(&setpnts, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Bcast failed");
  int step;
  step = setctrs % pool ? setctrs / pool + 1 : setctrs / pool;
  BcastI2d(meta, setctrs, 5);
  BcastD2d(data, setpnts, 4);
  if ( rank == pariah)
    if (stout > 1) PrintD2d(data, setpnts, 4); // complete data set 

  double **alldsts;                      // for gathering all distances
  CallocD2d(&alldsts, setctrs, setctrs); // note the '&'
  double **mydsts;                       // part of alldsts
  CallocD2d(&mydsts, step, setctrs);     // note the '&'

  ////////////////// PAIR-WISE CLOUD MATCHING /////////////////////////////////
  Spider(data, data, meta, meta, mydsts, wgtme, setctrs);
  /////////////////////////////////////////////////////////////////////////////

  // GATHERING //
  MPI_Gather(*mydsts, step * setctrs, MPI_DOUBLE, *alldsts, step * setctrs,
             MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // OUTPUT //
  if (rank == MASTER) {
    FILE *fprt;
    fprt = fopen(tabfn, "w");
    char *fmt = "%5d %5d %5d %5d %3d %3d %1.3lf %1.3lf\n";
    for (int i = 0; i < setctrs; i++)
      for (int j = 0; j < setctrs; j++) {
        double dst = alldsts[i][j];
        double tdst = alldsts[j][i];
        double dev = fabs(dst - tdst) / (dst + tdst);
        double mdst = (dst + tdst) / 2;
        fprintf(fprt, fmt,                                             //
                i, j, meta[i][1], meta[j][1], meta[i][0], meta[j][0], dst, dev);
      }
    fclose(fprt);
  }
  fflush(stdout);
  FreeD2d(&mydsts);
  FreeI2d(&meta);
  if (rank == MASTER)
    FreeD2d(&alldsts);
  if (stout) printf("%d:%s:%s:%s: bye,bye...\n", rank, __func__, __DATE__, __TIME__);
  MPI_Finalize();
  // avoid any code here - unpredictable behavior
  return 0;}
