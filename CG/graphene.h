// -----------------------------------------------------------------
// Header file for graphene on L^2 x Nt lattice
#ifndef GRAPHENE_H
#define GRAPHENE_H
#define USELAPACK

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <time.h>
#include <math.h>
#include "ranlxd.h"

#define r(...) creal(__VA_ARGS__)
#define i(...) cimag(__VA_ARGS__)
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Program parameters
typedef struct {
  int L;            // Length of each side of graphene lattice
  int t;            // Number of points in time direction
  int N;            // Total number of lattice sites
  double ka;        // Hopping parameter in action
  double beta;      // Inverse temperature
  double g;         // Coulomb strength
  double r0;        // Coulomb self-site effective radius
  double res;       // CG residual
} params_t;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// LAPACK routines
// http://www.netlib.org/lapack/explore-3.1.1-html/dgetrf.f.html
// http://www.netlib.org/lapack/explore-3.1.1-html/dgetri.f.html
void dgetrf_(int *N1, int *N2, double *store, int *lda, int *ipiv, int *stat);
void dgetri_(int *N, double *store, int *lda, int *ipiv, double *work,
             int *Nwork, int *stat);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void initR(int N, complex *R);
void CG(params_t info, double *phi, complex **D,
        complex *in, complex *out, int **hop);

#endif // GRAPHENE_H
// -----------------------------------------------------------------
