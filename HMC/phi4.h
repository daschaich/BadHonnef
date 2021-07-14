// -----------------------------------------------------------------
// Header file for phi^4 theory on L^D lattice
#ifndef PHI4_H
#define PHI4_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ranlxd.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Program parameters
typedef struct {
  int D;            // Number of dimensions
  int L;            // Length of lattice in each dimension
  int N;            // Total number of lattice sites, N = L^D
  double ka;        // Hopping parameter in action
  double la;        // Self-coupling in action
  int X;            // Use Omelyan integrator if nonzero
  double xi;        // Omelyan tunable parameter
  int ntraj;        // Number of trajectories
  double tlength;   // Trajectory length (not necessarily integer)
  int nstep;        // Number of steps per trajectory
} params_t;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int hmc(params_t info, double *phi, int **hop, double *mom);

double calcMag(int N, double *phi);
#endif // PHI4_H
// -----------------------------------------------------------------
