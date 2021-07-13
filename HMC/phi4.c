// -----------------------------------------------------------------
// Main code for phi^4 theory on L^D lattice
// Struct params_t defined in the header file phi.h
// Other functions implemented in utils.c
#include "phi4.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Usage info
void print_usage(char *argv0) {
  char *fmt = " %-9s %s\n";

  printf("\n%s [options]\n", argv0);
  printf("options:\n");
  printf(fmt, "D #", "dimension");
  printf(fmt, "L #", "length");
  printf(fmt, "k #", "kappa");
  printf(fmt, "l #", "lambda");
  printf(fmt, "X #", "use Omelyan integrator if nonzero");
  printf(fmt, "x #", "Omelyan tunable parameter");
  printf(fmt, "N #", "number of trajectories");
  printf(fmt, "t #", "trajectory length");
  printf(fmt, "s #", "number of steps per trajectory");
  printf(fmt, "r #", "random number seed");
  printf("\n");
  exit(1);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up neighbors for a D dimensional lattice of size N = L^D
// The index of (n_0, n_1, ..., n_{D - 1}) is i = sum_k n_k L^k
// The index of its neighbor in positive direction mu is hop[i][mu]
// The index of its neighbor in negative direction mu is hop[i][D + mu]
void initHop(params_t info, int **hop) {
  int D = info.D, L = info.L, N = info.N;
  int x, y, Lk;
  int xk, k, dxk;

  // Go through all the points
  for (x = 0; x < N; x++) {
    Lk = N;
    y  = x;

    // Go through the components k
    for (k = D - 1; k >= 0; k--) {
      Lk /= L;                        // pow(L, k)
      xk = y / Lk;                    // kth component
      y  = y - xk * Lk;               // y<-y%Lk

      // Forward
      if (xk < L - 1)
        dxk = Lk;
      else
        dxk = Lk * (1 - L);
      hop[x][k] = x + dxk;

      // Backward
      if (xk > 0)
        dxk = -Lk;
      else
        dxk = Lk * (L - 1);
      hop[x][k + D] = x + dxk;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char* argv[]) {
  params_t info;
  char c, *argp;
  int i, accept, seed;
  int **hop;
  double dtime, m = 0.0;
  double *phi, *mom;              // Field and conjugate momentum

  // Set defaults, to be overwritten by user
  info.D = 3;
  info.L = 6;
  info.ka = 0.185825;
  info.la = 1.1689;
  info.X = 1;
  info.xi = 0.1931833;
  info.ntraj = 10;
  info.tlength = 1.0;
  info.nstep = 20;
  seed = 42;

  // Load parameters from command line and put in struct for easy passing
  for (i = 1; i < argc; i++) {
    argp = &argv[i][1];
    c = argv[i][0];
    if (argv[i][1] == '\0') {
      if (i + 1 < argc) {
        argp = argv[i + 1];
        i++;
      }
      else
        argp = NULL;
    }
    switch(c) {
      case 'D' : info.D = atoi(argp); break;
      case 'L' : info.L = atoi(argp); break;
      case 'k' : info.ka = atof(argp); break;
      case 'l' : info.la = atof(argp); break;
      case 'X' : info.X = atoi(argp); break;
      case 'x' : if (info.X != 0) info.xi = atof(argp); break;
      case 'N' : info.ntraj = atoi(argp); break;
      case 't' : info.tlength = atof(argp); break;
      case 's' : info.nstep = atoi(argp); break;
      case 'r' : seed = atoi(argp); break;
      default : print_usage(argv[0]);
    }
  }

  // Total lattice volume
  info.N = info.L;
  for (i = 1; i < info.D; i++)
    info.N *= info.L;

  // Print parameters
  printf("D = %d\n", info.D);
  printf("L = %d\n", info.L);
  printf("N = %d\n", info.N);
  printf("kappa = %g\n", info.ka);
  printf("lambda = %g\n", info.la);
  printf("X = %d\n", info.X);
  if (info.X == 0)
    printf("Using Verlet integrator\n");
  else
    printf("Using Omelyan integrator with xi = %g\n", info.xi);
  printf("ntraj = %d\n", info.ntraj);
  printf("traj_length = %g\n", info.tlength);
  printf("nstep = %d\n", info.nstep);
  printf("eps = %g\n", info.tlength / (double)info.nstep);
  printf("seed = %d\n", seed);

  // Start timer
  dtime = -1.0 * (double)clock();

  // Allocate memory for fields, neighbors list, measurements
  phi = malloc(sizeof *phi * info.N);
  mom = malloc(sizeof *mom * info.N);
  hop = malloc(info.N * sizeof(int *));
  for (i = 0; i < info.N; i++)
    hop[i] = malloc(2 * info.D * sizeof(int));

  // Set up neighbors list 'hop'
  initHop(info, hop);

  // Initialize random number generator and phi field
  rlxd_init(1, seed);
  ranlxd(phi, info.N);           // Random initialization

  // Loop over trajectories
  printf("\nStarting trajectories\n");
  for (i = 0; i < info.ntraj; i++) {
    accept = hmc(info, phi, hop, mom);

    // Avoid unnecessary re-calculation of unchanged magnetization
    if (accept == 1 || i == 0)
      m = calcMag(info.N, phi);

    printf("MAG %.6g\n", m);
  }

  dtime += (double)clock();
  dtime /= CLOCKS_PER_SEC;
  printf("\n%d trajectories done in %.2g sec\n", info.ntraj, dtime);

  return 0;
}
// -----------------------------------------------------------------
