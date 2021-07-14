// -----------------------------------------------------------------
// Initialization, HMC and measurement routines for phi^4 theory
#include "phi4.h"
//#define REVERSE

static double twopi = 6.2831847307179586;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fill momentum field with normally distributed random numbers
void initMom(int N, double *mom) {
  // TODO: To be implemented

  // Here is a sample conversion from uniform to gaussian random numbers
#if 0
  double u[2], rootLog, g[2];
  ranlxd(u, 2);

  rootLog = sqrt(-2.0 * log(1.0 - u[0]));
  g[0] = rootLog * cos(twopi * u[1]);
  g[1] = rootLog * sin(twopi * u[1]);
  printf("%.6g, %.6g\n", g[0], g[1]);
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update the phi field by step epsilon
void updatePhi(int N, double eps, double *phi, double *mom) {
  // TODO: To be implemented
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculate force term and update the momentum field by step epsilon
void updateMom(params_t info, double eps, double *phi, int **hop,
               double *mom) {

  // TODO: To be implemented
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Single step of MD integrator
void step(params_t info, double *phi, int **hop, double *mom) {
  if (info.X == 0) {                    // Verlet
    // TODO: To be implemented
  }
  else {                                // Omelyan
    // TODO: To be implemented
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute the action, including the momentum contribution
//    S = Sum_x [ pi_x^2 / 2 + phi_x^2 + la * (phi_x^2 - 1)^2
//                - 2 ka * Sum_mu phi_x phi_{x + mu} ]
double calcAct(params_t info, double *phi, int **hop, double *mom) {
  // TODO: To be implemented
  return -99.0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute the change in the action
// Reduce roundoff errors by summing site-by-site differences
// (Shouldn't matter much for the small systems testable by this code)
double calcDelta(params_t info, double *phi, double *newPhi, int **hop,
                 double *mom, double *newMom) {

  // TODO: To be implemented
  return 99.0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Run a full HMC trajectory.
// Return 1 if accepted, 0 if rejected
int hmc(params_t info, double *phi, int **hop, double *mom) {
  int i, accept;
  double S, delta, rand;

  // Momentum heat bath
  initMom(info.N, mom);

  // Allocate and initialize copies to be updated
  double *newPhi = malloc(sizeof *newPhi * info.N);
  double *newMom = malloc(sizeof *newMom * info.N);
  for (i = 0; i < info.N; i++) {
    newPhi[i] = phi[i];
    newMom[i] = mom[i];
  }

  // Calculate initial action
  S = calcAct(info, phi, hop, mom);

  // Calculate MD trajectory
  for (i = 0; i < info.nstep; i++)
    step(info, newPhi, hop, newMom);

  // Accept / reject test
  delta = calcDelta(info, phi, newPhi, hop, mom, newMom);
  ranlxd(&rand, 1);
  if (delta < 0 || exp(-delta) > rand) {
    printf("ACCEPT: ");
    accept = 1;
    for (i = 0; i < info.N; i++)
      phi[i] = newPhi[i];
  }
  else {
    accept = 0;
    printf("REJECT: ");
  }
  printf("delta = %.6g start = %.6g end = %.6g\n", delta, S, S + delta);
  printf("EXP %.6g\n", exp(-delta));

#ifdef REVERSE
  // Optionally check reversibility
  // TODO: To be implemented
#endif

  free(newPhi);
  free(newMom);
  return accept;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute the magnetization density m = (Sum_x phi_x) / V
double calcMag(int N, double *phi) {
  // TODO: To be implemented
  return -99.0;
}
// -----------------------------------------------------------------
