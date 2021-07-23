// -----------------------------------------------------------------
// Initialization, HMC and measurement routines for phi^4 theory
#include "phi4.h"
//#define REVERSE

static double twopi = 6.2831847307179586;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fill momentum field with normally distributed random numbers
void initMom(int N, double *mom) {
  int i;
  double x, rootLog;

  ranlxd(mom, N);

  // Convert from uniform to gaussian random numbers, by pairs
  for (i = 0; i < N - 1; i += 2) {
    x = mom[i + 1];
    rootLog = sqrt(-2.0 * log(1.0 - mom[i]));
    mom[i] = rootLog * cos(twopi * x);
    mom[i + 1] = rootLog * sin(twopi * x);
  }

  // Handle case of odd number of sites
  if (N % 2 == 1) {
    ranlxd(&x, 1);
    rootLog = sqrt(-2.0 * log(1.0 - mom[N - 1]));
    mom[N - 1] = rootLog * cos(twopi * x);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update the phi field by step epsilon
void updatePhi(int N, double eps, double *phi, double *mom) {
  for (int i = 0; i < N; i++)
    phi[i] += eps * mom[i];
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculate force term and update the momentum field by step epsilon
void updateMom(params_t info, double eps, double *phi, int **hop,
               double *mom) {

  int i, mu;
  double phiCubed = 0, F = 0;

  for (i = 0; i < info.N; i++) {
    // Calculate force term for each site
    phiCubed = phi[i] * phi[i] * phi[i];
    F = 2.0 * phi[i] + 4.0 * info.la * (phiCubed - phi[i]);
    for (mu = 0; mu < info.D; mu++)
      F -= 2.0 * info.ka * (phi[hop[i][mu]] + phi[hop[i][info.D + mu]]);

    // Actually calculated -F above; subtract to fix sign
    mom[i] -= eps * F;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Single step of MD integrator
void step(params_t info, double *phi, int **hop, double *mom) {
  double eps = info.tlength / info.nstep;

  if (info.X == 0) {                    // Verlet
    updatePhi(info.N, 0.5 * eps, phi, mom);
    updateMom(info, eps, phi, hop, mom);
    updatePhi(info.N, 0.5 * eps, phi, mom);
  }
  else {                                // Omelyan
    updatePhi(info.N, info.xi * eps, phi, mom);
    updateMom(info, 0.5 * eps, phi, hop, mom);
    updatePhi(info.N, eps * (1.0 - 2.0 * info.xi), phi, mom);
    updateMom(info, 0.5 * eps, phi, hop, mom);
    updatePhi(info.N, info.xi * eps, phi, mom);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute the action, including the momentum contribution
//    S = Sum_x [ pi_x^2 / 2 + phi_x^2 + la * (phi_x^2 - 1)^2
//                - 2 ka * Sum_mu phi_x phi_{x + mu} ]
double calcAct(params_t info, double *phi, int **hop, double *mom) {
  int i, mu;
  double S = 0.0, phiSq, td;

  // Loop over all sites
  for (i = 0; i < info.N; i++) {
    // Momentum term
    S += 0.5 * mom[i] * mom[i];

    // Kinetic term, with td = Sum_mu phi_{x + mu}
    td = phi[hop[i][0]];
    for (mu = 1; mu < info.D; mu++)
      td += phi[hop[i][mu]];

    S -= 2.0 * info.ka * phi[i] * td;

    // Mass and self-interaction term
    phiSq = phi[i] * phi[i];
    S += phiSq + info.la * (phiSq - 1.0) * (phiSq - 1.0);
  }
  return S;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute the change in the action
// Reduce roundoff errors by summing site-by-site differences
// (Shouldn't matter much for the small systems testable by this code)
double calcDelta(params_t info, double *phi, double *newPhi, int **hop,
                 double *mom, double *newMom) {

  int i, mu;
  double delta = 0.0;
  double td, new_td, phiSq, newPhiSq;

  // Calculate the change in energy at each site
  for (i = 0; i < info.N; i++) {
    // Momentum term
    delta += 0.5 * (newMom[i] * newMom[i] - mom[i] * mom[i]);

    // Kinetic term
    td = phi[hop[i][0]];
    new_td = newPhi[hop[i][0]];
    for (mu = 1; mu < info.D; mu++) {
      td += phi[hop[i][mu]];
      new_td += newPhi[hop[i][mu]];
    }
    delta -= 2.0 * info.ka * (newPhi[i] * new_td - phi[i] * td);

    // Mass and self-interaction term
    phiSq = phi[i] * phi[i];
    newPhiSq = newPhi[i] * newPhi[i];
    td = (newPhiSq - 1.0) * (newPhiSq - 1.0) - (phiSq - 1.0) * (phiSq - 1.0);
    delta += newPhiSq - phiSq + info.la * td;
  }
  return delta;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Run a full HMC trajectory.
// Return 1 if accepted, 0 if rejected
int hmc(params_t info, double *phi, int **hop, double *mom) {
  int i, accept;
  double S, delta, rand, td;

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
  td = exp(-delta);
  ranlxd(&rand, 1);
  if (delta < 0 || td > rand) {
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
  printf("EXP %.6g\n", td);

#ifdef REVERSE
  // Optionally check reversibility
  for (i = 0; i < info.N; i++)
    newMom[i] = -1.0 * newMom[i];

  for (i = 0; i < info.nstep; i++)
    step(info, newPhi, hop, newMom);

  // For simple test, take difference of total action
  td = calcAct(info, newPhi, hop, newMom);
  printf("REVERSED: delta = %.16g start = %.16g end = %.16g\n",
         td - S, S, td);
#endif

  free(newPhi);
  free(newMom);
  return accept;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute the magnetization density m = (Sum_x phi_x) / V
double calcMag(int N, double *phi) {
  int i;
  double M = 0.0;

  // Loop over all sites
  for (i = 0; i < N; i++)
    M += phi[i];

  return (M / (double)N);
}
// -----------------------------------------------------------------
