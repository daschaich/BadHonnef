// -----------------------------------------------------------------
// Includes initialization functions for simulations of graphene.
#include "graphene.h"

static double twopi = 6.2831847307179586;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fill complex vector with normally distributed random numbers
void initR(int N, complex *R) {
  int i;
  double u[2], rootLog;

  // Convert from uniform to gaussian random numbers, by pairs
  for (i = 0; i < N; i++) {
    ranlxd(u, 2);
    rootLog = sqrt(-2.0 * log(1.0 - u[0]));
    R[i] = rootLog * cos(twopi * u[1]);        // Real part
    R[i] += I * rootLog * sin(twopi * u[1]);   // Imaginary part
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Solve chi = (D.Ddag)^{-1} (PF)
void CG(params_t info, double *phi, complex **D, complex *in,
        complex *out, int **hop) {

  int i, j, k, N = info.N, m = 2 * info.L * info.L;
  double alpha, beta, denom, gbt = info.beta * info.g / info.t;
  double rsq = 0.0, newRsq = 0.0;
  complex tc, *res, *vec, **DDdag;

  // Allocations and initializations
  res = malloc(sizeof *res * N);
  vec = malloc(sizeof *vec * N);
  DDdag = malloc(N * sizeof(complex *));
  for (i = 0; i < N; i++) {
    DDdag[i] = malloc(N * sizeof(complex));
    res[i] = in[i];
    vec[i] = in[i];
    out[i] = 0;
    j = hop[i][4];
    D[i][j] -= I * gbt * phi[j] * (2 * floor((i - m + N) / N) - 1);
  }

  // Construct DDdag
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      DDdag[i][j] = 0;
      for (k = 0; k < N; k++)
        DDdag[i][j] += D[i][k] * conj(D[j][k]);
    }
  }

  // Will exit on convergence or warn if all N iterations exhausted
  for (k = 0; k < N; k++) {
    // Update out and get new residual using current residual
    rsq = 0;
    denom = 0;
    for (i = 0; i < N; i++) {
      rsq += res[i] * conj(res[i]);
      for (j = 0; j < N; j++)
        denom += conj(vec[i]) * DDdag[i][j] * vec[j];
    }
    alpha = rsq / denom;
    for (i = 0; i < N; i++) {
      out[i] += alpha * vec[i];
      tc = DDdag[i][0] * vec[0];
      for (j = 1; j < N; j++)
        tc += DDdag[i][j] * vec[j];

      res[i] -= alpha * tc;
    }

    // Exit if new residual is small enough
    newRsq = res[0] * conj(res[0]);;
    for (i = 1; i < N; i++)
      newRsq += res[i] * conj(res[i]);
    if (newRsq < info.res * info.res) {
      printf("CG converged with rsq = %g ", newRsq);
      printf("after %d iterations\n", k + 1);
      break;
    }

    // Update vec using new residual
    beta = newRsq / rsq;
    printf("beta%d: %g = %g / %g\n", k + 1, beta, newRsq, rsq);
    for (i = 0; i < N; i++) {
      vec[i] *= beta;
      vec[i] += res[i];
    }
  }

  // Warn if Krylov subspace should cover full dimensionality
  if (k == N) {
    printf("WARNING: %d CG iterations exhausted ", N);
    printf("with rsq = %g\n", newRsq);
  }

  // Clean up
  for (i = 0; i < N; i++) {
    j = hop[i][4];
    D[i][j] += I * gbt * phi[j] * (2 * floor((i - m + N) / N) - 1);
    free(DDdag[i]);
  }
  free(vec);
  free(res);
  free(DDdag);
}
// -----------------------------------------------------------------
