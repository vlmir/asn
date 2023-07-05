
#include "icp.h"
/// for both icp.c asine.c
// vars below cause linker errors
const int D = 3; // dimentionality, xyz in this case; TODO: as arg
const int MAXITERS = 512;

/// only for icp.c
const double MAXDBL = 1E9;        // for initializing dms vars
const double TAUSQ = DBL_EPSILON; // 2^-53 ~~ 2.22*e-16
// const double TAU = 1.49e-8;  // tolerance for ICP computation; ori:1e-9;
// int x2z[64];
///////////////////////////////////////////////////////////////////////////////

int Twist(double RR[24][D][D], double CVX[D][D], double CVY[D][D], Ret *auxp) {
  // Attn: no closest points here !!
  /*
   * computes 24 det-positive rotation matrices R for 2 co-var matrices
   * R=P*D*Vz*(Vw)^T ~ Va*(Vb)^T :: CahillND2000
   * 6 P[3][3] times 8 D[3][3] / 2
   */
  /*
  * 48 matreces spanning all possible rotations and reflections
  * 24 with det == +1, 24 with det == -1
  */
  // 6 3*3 permutation matrices
  const int PERMUT[18] = {// the 1st is I3, the 4th its transpose
                          0, 1, 2, 2, 0, 1, 1, 2, 0,  // det == +1
                          2, 1, 0, 1, 0, 2, 0, 2, 1}; // det == -1
  // 8 3*3 diagonal matrices
  const int DIAGON[24] = {
      1, 1, 1,  1, -1, -1, -1, 1, -1, -1, -1, 1, // det == +1
      1, 1, -1, 1, -1, 1,  -1, 1, 1,  -1, -1, -1 // det == -1
  };

  /// co-variance component ///
  double AV[D][D]; // eigen vectors for nefosX; ori:evec
  double BV[D][D]; // eigen vectors for nefosY; ori:evecp
  double dummy[D];
  auxp->twist[0] += jacobi(D, CVX, dummy, AV, MAXITERS);
  auxp->twist[1] += jacobi(D, CVY, dummy, BV, MAXITERS);
  double TV[D][D];   // transpose of AV
  double BVTV[D][D]; // null matrix if A==B ?
  Transpose(D, D, AV, D, D, TV);
  mulmat3(D, D, BV, D, D, TV, D, D, BVTV); // TODO firs matrix TV <-> BV ?
  double det = Determinant(D, BVTV);
  if (det == 0)
    return 0;
  double PM[D][D];

  /// modifier matrices P*D ///
  int ini, fin; // DIAGON indeces
  int r = 0;
  for (int i = 0; i < 6; i++) { // 6 permutation matrices
    // constructing one of the 6 permutation matreces PM
    int n = i * D; // 0..15, offset
    for (int j = 0; j < D; j++) {
      for (int k = 0; k < D; k++) {
        if (k == PERMUT[n])
          PM[j][k] = 1;
        else
          PM[j][k] = 0;
      }
      n++;
    }
    // double dpm = Determinant(D, PM);
    // printf("r:%d i:%d dpm:%lf \n", r, i, dpm);

    // for each permutatin matrix PM eight PD products
    double PD[D][D];
    NullDm(D, D, PD);

    if ((det > 0 && i < 3) || (det < 0 && i > 2)) { // det PM > 0
      ini = 0;
      fin = 12;
    } else if ((det > 0 && i > 2) || (det < 0 && i < 3)) { // det PM < 0
      ini = 12;
      fin = 24;
    }
    for (int m = ini; m < fin; m += 3) { // 4 D matrices
      for (int k = 0; k < D; k++)
        for (int l = 0; l < D; l++)
          PD[k][l] = PM[k][l] * DIAGON[m + k];

      // double dpd = Determinant(D, PD);
      // printf("r:%d m:%d dpd:%lf \n", r, m, dpd);
      // PrintDm(D, D, PD);

      double R[D][D];
      mulmat3(D, D, PD, D, D, BVTV, D, D, R);
      double dr = Determinant(D, R);
      if (dr < 0)
        fprintf(stderr, "Negative det R: r:%d dr:%lf \n", r, dr);
      // segregating rotation and reflection matrices ()
      // copying only det positive matrices into RR
      for (int k = 0; k < D; k++)
        for (int l = 0; l < D; l++)
          RR[r][k][l] = R[k][l];
      r++;
    }
  } // i<6
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

/*
       * ICP point matching
 * x, y, z: clouds of d-dimentional points
       * returns array of shortest squared dsts for cloud y in cloud x
 * m: size of x
 * n: size of y and z
 */
double *OptimeProximo(int d, int m, int n, double **x, double **y, double **z,
                      Ret *auxp) {
  for (int i = 0; i < m; i++)
    auxp->optimeproximo[i] = 0;
  double *y2z = NULL; // sum of squares between y and z points
  if ((y2z = (double *)calloc(n, sizeof(double))) == NULL)
    Error("calloc failed");
  for (int j = 0; j < n; j++) {
    double minss; // ori:min_dist^2
    minss = MAXDBL;
    int ij = -1; // index corresponding to mindist
    for (int i = 0; i < m; i++) {
      double ss = EDsquare(d, y[j], x[i]);
      if (ss < minss) {
        minss = ss;
        ij = i;
      }
    }
    for (int k = 0; k < d; k++)
      z[j][k] = x[ij][k];
    y2z[j] = minss;
    auxp->optimeproximo[ij]++; // count how many times each x pnt is used
  }
  return y2z; // unsorted
}

///////////////////////////////////////////////////////////////////////////////

/// ICP point registration ///
// nuageX, nuageY: point clouds, mass centered -> no translations
double PuntaOptima(double **nuageX, double **nuageY, int ctrXpnts, int ctrYpnts,
                   int wgtme, Ret *auxp) {
  const int Q = 4;        // quaternion size
  double mseMin = MAXDBL; // returned value sqrt(mseMin)
  double mseNee = MAXDBL; // mseMin backup
  double delta = MAXDBL;  // ALWAYS positive; ori:dmsdiff
  const double TAU = sqrt(TAUSQ);
  double mse; // should stay here; TODO why?? // weighted if wgtme!=0

  /// heap vars ///
  double *src2cp; // sqared dsts between src and closest points in trg,
                  // ori:min2_dist^2
  if ((src2cp = (double *)calloc(ctrYpnts, sizeof(double))) == NULL)
    Error("calloc failed");

  // PuntaOptima//---------------------- main loop
  // ------------------------------------///
  // auxp->puntaoptima[0] = 0; // number of closest point iterations
  int cpis; // number of closest point iterations
  int jrotsQ = auxp->puntaoptima[1];
  for (cpis = 0; cpis < MAXITERS; cpis++) { // normally never reached
    // ori condition:
    // if((dmsdiff<tol && dmsdiff>=0) || cpiss>=500){ break; }
    double **nuageZ; // target points closest to sourc points, ori: cp
    // used for XV computation at the start of the main loop
    CallocD2d(&nuageZ, ctrYpnts, D);
    double **nuageR;                 // used only in the end of main loop;
    CallocD2d(&nuageR, ctrYpnts, D); // nuageR = nuageY*QT

    /// 1. CLOSEST POINTS ///

    /* a.k.a.:
     * data association
     * point matching
     * correspondence finding
     */
    // building nuageZ and src2cp (src2cp computed also in MSerror TODO: decide)
    // ctrXpnts used only here
    src2cp = OptimeProximo(D, ctrXpnts, ctrYpnts, nuageX, nuageY, nuageZ, auxp);
    CentrifyD2d(nuageZ, ctrXpnts, D);

    /// 2. REGISTRATION ///
    // computing weighted mse//
    // each ctr must have at least 2 ipnts
    mse = MSerror(D, ctrYpnts, nuageY, nuageZ, wgtme);
    if (mse < mseMin)
      mseMin = mse;
    ///////////////////////////////////////////////////////////////////////////
    // TODO double check the logics !!!
    // the only assign to delta; the only use of mseNee
    delta = fabs(mseNee - mse);
    if ((delta < TAU)
        //&& (mse <= mseMin)//PuntaOptima: reached the limit of iterations:500
        // delta:0.000000
        )
      break;      // !!!
    mseNee = mse; // for the next iteration, the only assign

    ///////////////////////////////////////////////////////////////////////////

    /// Registration ///

    //<
    double XV[D][D]; // cross-covariance matrix, sum of 3d mtxs; ori:Spx
    double trace;    // ori:TrA
    double qr[Q];    // rotation quaternion, corresponds to max eigen value
    double XT[D][D]; // transpose of XV, ori: tatr
    double QX[Q][Q]; // quaternion matrix; ori:Q
    double qv[Q];    // eigenvalues of quaternion matrix; ori:val
    double QV[Q][Q]; // eigenvectors of quaternion matrix; ori:vec
    double QR[D][D]; // ratation matrix; ori:ROT
    double QT[D][D]; // QR transposed; ori:MTM
    //> no change if the vars between <> declared in or out the main loop

    /* not needed with mass-centered data
double meanA[D], meanB[D], meanC[D], mean0[D];
double qt[D];   // translation vector, no need with mass-centered data
     */

    /// computing XV ///
    // Note: no difference nuageY-nuageZ or nuageZ-nuageY, sym matrix
    Xwash(ctrYpnts, D, nuageY, nuageZ, XV);
    /// make quaternion matrix ///
    Transpose(D, D, XV, D, D, XT);
    trace = XV[0][0] + XV[1][1] + XV[2][2];
    QX[0][0] = trace;
    QX[1][1] = XV[0][0] + XT[0][0] - trace;
    QX[2][2] = XV[1][1] + XT[1][1] - trace;
    QX[3][3] = XV[2][2] + XT[2][2] - trace;
    QX[0][1] = XV[1][2] - XT[1][2];
    QX[0][2] = XV[2][0] - XT[2][0];
    QX[0][3] = XV[0][1] - XT[0][1];
    QX[1][0] = XV[1][2] - XT[1][2];
    QX[2][0] = XV[2][0] - XT[2][0];
    QX[3][0] = XV[0][1] - XT[0][1];
    QX[1][2] = XV[0][1] + XT[0][1];
    QX[1][3] = XV[0][2] + XT[0][2];
    QX[2][1] = XV[1][0] + XT[1][0];
    QX[2][3] = XV[1][2] + XT[1][2];
    QX[3][1] = XV[2][0] + XT[2][0];
    QX[3][2] = XV[2][1] + XT[2][1];

    // find eigen values and vectors
    jrotsQ += jacobi(Q, QX, qv, QV, MAXITERS);

    // largest eigen value of quaternion matrix //
    NullD1d(Q, qr); // ini to nils
    qr[0] = 1.;
    double maxqval = 0.;
    int maxi = 0;
    for (int i = 0; i < Q; i++) {
      qr[i] = 0.;
      if (qv[i] > maxqval) {
        maxqval = qv[i];
        maxi = i;
      }
    }
    // quaternion corresponds to the maximal eigen value of quaternion matrix
    for (int i = 0; i < Q; i++)
      qr[i] = QV[i][maxi];

    // define rotation matrix
    QR[0][0] = qr[0] * qr[0] + qr[1] * qr[1] - qr[2] * qr[2] - qr[3] * qr[3];
    QR[1][1] = qr[0] * qr[0] - qr[1] * qr[1] + qr[2] * qr[2] - qr[3] * qr[3];
    QR[2][2] = qr[0] * qr[0] - qr[1] * qr[1] - qr[2] * qr[2] + qr[3] * qr[3];
    QR[0][1] = 2. * (qr[1] * qr[2] - qr[0] * qr[3]);
    QR[0][2] = 2. * (qr[1] * qr[3] + qr[0] * qr[2]);
    QR[1][0] = 2. * (qr[1] * qr[2] + qr[0] * qr[3]);
    QR[1][2] = 2. * (qr[2] * qr[3] - qr[0] * qr[1]);
    QR[2][0] = 2. * (qr[1] * qr[3] - qr[0] * qr[2]);
    QR[2][1] = 2. * (qr[2] * qr[3] + qr[0] * qr[1]);

    /// 3. APPLYING REGISTRATION ///
    // define nuageR nuageR = nuageY*QT
    Transpose(D, D, QR, D, D, QT);
    // NullD2d(ctrYpnts, D, nuageR); // if nuageR allocated outside the loop
    for (int b = 0; b < ctrYpnts; b++)
      for (int l = 0; l < D; l++)
        for (int k = 0; k < D; k++)
          nuageR[b][l] += nuageY[b][k] * QT[k][l];
    // update nuageY
    for (int b = 0; b < ctrYpnts; b++)
      for (int l = 0; l < D; l++)
        nuageY[b][l] = nuageR[b][l]; // nuageY<-nuageR
    cpis++;
    FreeD2d(&nuageZ);
    FreeD2d(&nuageR);
  }
  /// ------------------ end of main loop --------------------------------///
  auxp->puntaoptima[0] = cpis;
  auxp->puntaoptima[1] = jrotsQ;

  if (cpis >= MAXITERS)
    fprintf(stderr, //
            "PuntaOptima: reached the limit of iterations:%d delta:%lf\n", cpis,
            delta);
  free(src2cp);
  return sqrt(mseMin);
} // end of PuntaOptima()

/////////////////////////////////////////////////////////////////////////

/*
 * computes tentative GLOBAL registration for 2 point clouds
 * with 24 evenly spaced initial rotation matrices
 */
double RotMall(double **nefosX, double **nefosY, int ctrXpnts, int ctrYpnts,
               int wgtme, Ret *auxp) {
  // TODO supply func pointer as a parameter for the algoritm
  const double TAU = sqrt(TAUSQ);
  int test_rot = 0;
  double dstmin = MAXDBL; // the minimum over all R matrices tested; returened
  double RR[24][D][D];    // 24-available inital rotation matries; ori:R;
  double CVX[D][D];       // used only in Twist() for jacobi()
  Squash(ctrXpnts, D, nefosX, CVX);
  double CVY[D][D]; // used only in Twist() for jacobi()
  Squash(ctrYpnts, D, nefosY, CVY);
  Twist(RR, CVX, CVY, auxp);
  /////////////////////////////////////////////////////////////////////////
  // walking over the 24 possible rotation matrices
  // not known in advance which one will work best
  for (int r = 0; r < 24; r++) { // r == ori:iter
    if (dstmin < TAU)
      break;        // makes no sense to optimize further
    double R[D][D]; // immutable; ori:R?
    /// copying RR into R ///
    for (int i = 0; i < D; i++)
      for (int j = 0; j < D; j++) {
        R[i][j] = RR[r][i][j]; // RR and r used only here
      }

    if (test_rot) {
      double T[D][D];
      Transpose(D, D, R, D, D, T);
      double TR[D][D];
      mulmat3(D, D, T, D, D, R, D, D, TR);
      for (int k = 0; k < D; k++)
        if (fabs(TR[k][k] - 1.) > TAU)
          fprintf(stderr, "TR != I3: r:%d TR:%lf \n", r, TR[k][k]);
    }

    double **nefosR;                 // nefosY rotated, changes iteratively
    CallocD2d(&nefosR, ctrYpnts, D); // must be nullified
    // ATTN: aborted if nefosR defined outside of the R loop!!

    // nefosR = nefosY * R, pure rotation //
    for (int v = 0; v < ctrYpnts; v++)
      for (int j = 0; j < D; j++)
        for (int k = 0; k < D; k++)
          nefosR[v][j] += nefosY[v][k] * R[k][j];
    double dst = PuntaOptima(nefosX, nefosR, ctrXpnts, ctrYpnts, wgtme, auxp);

    if (dst < dstmin) {
      dstmin = dst;
      auxp->rotmall[0] = r + 1;
    }
    FreeD2d(&nefosR);
  } // r < 24
  return dstmin;
} // RotMall

/*
 * JACOBI.C  module
 * jacobi() originates from 'Numerical Recipes in C'
 * rotate() replaces the original macro ROTATE()
*/

typedef unsigned dimension;
typedef unsigned iterations;
/* rotate<() replaces the original macro ROTATE */
void rotate(int n, double S[n][n], int i, int j, int k, int l, double s,
            double tau) {
  double g = S[i][j];
  double h = S[k][l];
  S[i][j] = g - s * (h + g * tau);
  S[k][l] = h + s * (g - h * tau);
}

/* S: symmetric matrix
   n: dimensionality of S, MUST precede S in the argument list
   w: contains eigenvalues sorted in descending order
   V: contains coloumn eigenvectors of S properly sorted
   nrot: maximum number of rotations allowed in jacobi()
   RETURNS: >0 : number of jacobi iterations
            -1 : number of jacobi iterations exceeded MAXITERS
*/
int jacobi(int n, double S[n][n], double *w, double V[n][n],
           const int MAXITERS) {
  // iterations i, j, k, ip, iq; // ori: unsigned int
  int i, j, k, ip, iq;
  double thrash, theta, tau, t, sm, s, h, g, c;
  double p;
  double b[n];
  double z[n];
  int nrot;

  for (ip = 0; ip <= n - 1; ip++) {
    for (iq = 0; iq <= n - 1; iq++)
      V[ip][iq] = 0.0;
    V[ip][ip] = 1.0;
  }
  for (ip = 0; ip <= n - 1; ip++) {
    b[ip] = w[ip] = S[ip][ip];
    z[ip] = 0.0;
  }
  nrot = 0;

  for (i = 1; i <= MAXITERS; i++) {
    sm = 0.0;
    for (ip = 0; ip <= n - 2; ip++) {
      for (iq = ip + 1; iq <= n - 1; iq++)
        sm += fabs(S[ip][iq]);
    }
    if (sm == 0.0) {
      // eigenvalues & eigenvectors sorting
      for (i = 0; i < n - 1; i++) {
        p = w[k = i];
        for (j = i + 1; j <= n - 1; j++)
          if (w[j] >= p)
            p = w[k = j];
        if (k != i) {
          w[k] = w[i];
          w[i] = p;
          for (j = 0; j <= n - 1; j++) {
            p = V[j][i];
            V[j][i] = V[j][k];
            V[j][k] = p;
          }
        }
      }
      // restore symmetric matrix S
      for (i = 1; i <= n - 1; i++) {
        for (j = 0; j < i; j++)
          S[j][i] = S[i][j];
      }
      return (nrot);
    } // end of if sm
    thrash = i < 4 ? 0.2 * sm / (n * n) : 0.0;
    for (ip = 0; ip <= n - 2; ip++) {
      for (iq = ip + 1; iq <= n - 1; iq++) {
        g = 100.0 * fabs(S[ip][iq]);
        if (i > 4 && (fabs(w[ip]) + g) == fabs(w[ip]) &&
            (fabs(w[iq]) + g) == fabs(w[iq]))
          S[ip][iq] = 0.0;
        else if (fabs(S[ip][iq]) > thrash) {
          h = w[iq] - w[ip];
          if ((fabs(h) + g) == fabs(h))
            t = (S[ip][iq]) / h;
          else {
            theta = 0.5 * h / (S[ip][iq]);
            t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0)
              t = -t;
          }
          c = 1.0 / sqrt(1 + t * t);
          s = t * c;
          tau = s / (1.0 + c);
          h = t * S[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          w[ip] -= h;
          w[iq] += h;
          S[ip][iq] = 0.0;
          for (j = 0; j < ip; j++)
            rotate(n, S, j, ip, j, iq, s, tau);
          for (j = ip + 1; j < iq; j++)
            rotate(n, S, ip, j, j, iq, s, tau);
          for (j = iq + 1; j <= n - 1; j++)
            rotate(n, S, ip, j, iq, j, s, tau);
          for (j = 0; j <= n - 1; j++)
            rotate(n, V, j, ip, j, iq, s, tau);
          ++nrot;
        }
      }
    }
    for (ip = 0; ip <= n - 1; ip++) {
      b[ip] += z[ip];
      w[ip] = b[ip];
      z[ip] = 0.0;
    }
  }            // end of main for
  return (-1); /* Too many iterations in jacobi() */
} /* End of jacobi() */
