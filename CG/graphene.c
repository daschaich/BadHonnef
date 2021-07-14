// -----------------------------------------------------------------
// Main code for graphene on L^2 x Nt lattice
// Struct params_t defined in the header file graphene.h
// Other functions implemented in utils.c
#include "graphene.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Usage info
void print_usage(char *argv0) {
  char *fmt = " %-9s %s\n";

  printf("\n%s [options]\n", argv0);
  printf("options:\n");
  printf(fmt, "L #", "length of each side of graphene lattice");
  printf(fmt, "t #", "time");
  printf(fmt, "B #", "beta");
  printf(fmt, "k #", "kappa");
  printf(fmt, "g #", "Coulomb strength");
  printf(fmt, "R #", "Coulomb self-interaction effective radius");
  printf(fmt, "r #", "CG residual");
  printf(fmt, "S #", "random number seed");
  printf("\n");
  exit(1);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up neighbors for honeycomb lattice
// The first three components of each hop[x] contain the indices
//              of the point's three neighbors on the same timeslice
// The next two contain its forward and backward temporal neighbors
// The sixth component records whether point i is red or black
void initHop(params_t info, int **hop) {
  int L = info.L, N = info.N;
  int W = 2 * L;    // "Width"
  int m = W * L;    // Sites on each timeslice
  int i, x;

  // Go through all the points
  for (x = 0; x < N; x++) {
    i = x % m;   // Position within timeslice

    // "Previous" neighbor
    if (i % W == 0)               // Start of each "row"
      hop[x][0] = x + W - 1;
    else
      hop[x][0] = x - 1;

    // "Next" neighbor
    if ((i + 1) % W == 0)   // End of each "row"
      hop[x][1] = x + 1 - W;
    else
      hop[x][1] = x + 1;

    // Now for the fun part -- inter-"row" neighbors
    int row = (int)(floor(i / W)) % 2;

    if (row == 0) {       // Generic even "row"
      if (i % 2 == 0) {
        hop[x][2] = x + W;
        hop[x][5] = 0;    // "Red"
      }
      else {
        hop[x][2] = x - W;
        hop[x][5] = 1;    // "Black"
      }
    }
    else {                // Generic odd "row"
      if (i % 2 == 0) {
        hop[x][2] = x - W;
        hop[x][5] = 1;    // "Black"
      }
      else {
        hop[x][2] = x + W;
        hop[x][5] = 0;    // "Red"
      }
    }

    // Correct for boundary conditions
    // First row
    if (i < L && i % 2 == 1)
      hop[x][2] = x + m - L;
    else if (i < W && i % 2 == 1)
      hop[x][2] = x + m - W - L;
    // Last row depends on whether L is even or odd
    else if (i > (m - 1 - L) && L % 2 == 0 && i % 2 == 1)
      hop[x][2] = x + L - m;
    else if (i > (m - 1 - L) && L % 2 == 1 && i % 2 == 0)
      hop[x][2] = x + L - m;
    else if (i > (m - 1 - W) && L % 2 == 0 && i % 2 == 1)
      hop[x][2] = x + L + W - m;
    else if (i > (m - 1 - W) && L % 2 == 1 && i % 2 == 0)
      hop[x][2] = x + L + W - m;

    // Finally, forward and backward temporal neighbors
    if (x > (N - 1) - m)
      hop[x][3] = i;
    else
      hop[x][3] = x + m;

    if (x < m)
      hop[x][4] = N - m + i;
    else
      hop[x][4] = x - m;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Initialize the potential between two points in the unit cell
// Don't include image charges for now
void initV(params_t info, double **V, int **hop) {
  int i, j, k, l;
  int L = info.L;
  int W = 2 * L;    // "Width"
  int m = W * L;    // Sites on each timeslice
  int M = 1000;     // Temporary number of image charges in each direction
                    // Will eventually be chosen by convergence
  double invdelta = info.t / info.beta; // delta always in denominator
  double *x, *y, dx, dy, td, td2;
  double sqrt3ov2 = 0.5 * sqrt(3.0);
  double R2x = 1.5 * (double)L;
  double R2y = sqrt3ov2 * (double)L;
  double R1y = 2.0 * R2y;

  x = malloc(sizeof *x * m);
  y = malloc(sizeof *y * m);

  for (i = 0; i < m; i++) {
    x[i] = (double)(floor(i / W)) * 1.5  - 0.5 * (double)(hop[i][5]);
    y[i] = (i % W) * sqrt3ov2;
    if (y[i] >= 2.0 * sqrt3ov2 * (-x[i] + 2 * L)) {
      x[i] -= R2x;
      y[i] -= R2y;
    }
    if (y[i] <= 2.0 * sqrt3ov2 * (x[i] - L)) {
      x[i] -= R2x;
      y[i] += R1y - R2y;
    }
  }
  for (i = 0; i < m; i++) {
    V[i][i] = 2.0 * invdelta / info.r0;
    for (j = 0; j < i; j++) {
      dx = x[i] - x[j];
      dy = y[i] - y[j];
      V[i][j] = invdelta / sqrt(pow(dx, 2) + pow(dy, 2));
      for (k = -M; k < M + 1; k++) {
        for (l = -M; l < M + 1; l++) {
          if ((k != 0 || l != 0) && abs(k + j) <= M) {
            // Take out the zero mode
            // For the diagonal terms this cancels the image charges
            td = sqrt(pow(dx + l * R2x, 2) + pow(dy + k * R1y + l * R2y, 2));
            td2 = sqrt(pow(l * R2x, 2) + pow(k * R1y + l * R2y, 2));
            V[i][j] += invdelta * (1 / td - 1 / td2);
          }
        }
      }
    V[j][i] = V[i][j];
    }
  }
  free(x);
  free(y);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Initialize the D matrix without the terms containing phi
// As most values are zero, should eventually do away with this matrix
void initD(params_t info, complex **D, int **hop) {
  int i, j, N = info.N, L = info.L;
  int W = 2 * L;    // "Width"
  int m = W * L;    // Sites on each timeslice
  double g = info.g, ka = info.ka;
  double delta = info.beta / info.t;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
      D[i][j] = 0;

  }
  for (i = 0; i < N; i++) {
    D[hop[i][3]][hop[i][0]] = ka * delta * (1 - 2 * floor((i + m) / N));
    D[hop[i][3]][hop[i][1]] = D[hop[i][3]][hop[i][0]];
    D[hop[i][3]][hop[i][2]] = D[hop[i][3]][hop[i][0]];

    D[i][hop[i][4]] = (-1 + g * g * delta / info.r0)
                    * (2 * floor((i - m + N) / N) - 1);
    D[i][i] = 1;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Invert the potential V using LAPACK
// If LAPACK is not available, can only directly invert
//                         the 2 x 2 matrix for the two-site problem
void calcVInv(params_t info, double **VInv) {
  int i, j, m = 2 * info.L * info.L;
#ifdef USELAPACK
  int stat = 0;
  int *ipiv = malloc(sizeof *ipiv * m);   // Pivot indices
  double *store = malloc(sizeof *store * m * m);
  double *work = malloc(sizeof *work * m);

  printf("Using LAPACK to invert V\n");
  // Convert VInv to column-major form used by LAPACK
  for (i = 0; i < m; i++) {
    for (j = 0; j < m; j++)
      store[j * m + i] = VInv[i][j];
  }

  dgetrf_(&m, &m, store, &m, ipiv, &stat);
  dgetri_(&m, store, &m, ipiv, work, &m, &stat);

  // Copy back into row-major form
  for (i = 0; i < m; i++) {
    for (j = 0; j < m; j++)
      VInv[i][j] = store[j * m + i];
  }
  free(ipiv);
  free(store);
  free(work);
#else
  double V[2][2], det;

  // Only for two-site problem
  if (m != 2) {
    fprintf(stderr, "Error: need LAPACK for more than two sites\n");
    fflush(stderr);
    exit(1);
  }

  for (i = 0; i < m; i++) {
    for (j = 0; j < m; j++)
      V[i][j] = VInv[i][j];
  }

  det = V[0][0] * V[1][1] - V[0][1] * V[1][0];
  VInv[0][0] = V[1][1] / det;
  VInv[1][1] = V[0][0] / det;
  VInv[1][0] = -V[1][0] / det;
  VInv[0][1] = -V[0][1] / det;
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char* argv[]) {
  params_t info;
  char c, *argp;
  int i, j, seed, m, N;
  int **hop;
  double dtime, gbt;
  double *phi, *mom;              // Field and conjugate momentum
  double **VInv;                  // Coulomb potential inverse
  complex *R, *PF, *chi;          // Pseudofermions and related vectors
  complex **D;                    // Fermion operator (inefficient)

  // Set default parameters, to be overwritten by user
  info.L = 2;
  info.t = 8;
  info.beta = 16.0;
  info.ka = 0.25;
  info.g = 0.5;
  info.r0 = 0.1;
  info.res = 1e-8;
  seed = 42;

  // Load parameters from command line into struct for easy passing
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
      case 'L' : info.L = atoi(argp); break;
      case 't' : info.t = atoi(argp); break;
      case 'B' : info.beta = atoi(argp); break;
      case 'k' : info.ka = atof(argp); break;
      case 'g' : info.g = atof(argp); break;
      case 'R' : info.r0 = atof(argp); break;
      case 'r' : info.res = atof(argp); break;
      case 'S' : seed = atoi(argp); break;
      default : print_usage(argv[0]);
    }
  }

  m = 2 * info.L * info.L;  // Number of sites on each timeslice
  info.N = m * info.t;          // Total number of sites

  // Print parameters
  printf("L = %d\n", info.L);
  printf("t = %d\n", info.t);
  printf("N = %d\n", info.N);
  printf("beta = %g\n", info.beta);
  printf("kappa = %g\n", info.ka);
  printf("g = %g\n", info.g);
  printf("r0 = %g\n", info.r0);
  printf("res = %g\n", info.res);
  printf("seed = %d\n", seed);

  // Start timer
  dtime = -1.0 * (double)clock();

  // Allocate memory for fields and neighbors list
  phi = malloc(sizeof *phi * info.N);
  mom = malloc(sizeof *mom * info.N);
  hop = malloc(info.N * sizeof(int *));
  for (i = 0; i < info.N; i++)
    hop[i] = malloc(6 * sizeof(int));

  // Allocate matrices
  VInv = malloc(m * sizeof(double *));
  for (i = 0; i < m; i++)
    VInv[i] = malloc(m * sizeof(double));

  D = malloc(info.N * sizeof(complex *));
  for (i = 0; i < info.N; i++)
    D[i] = malloc(info.N * sizeof(complex));

  // Allocate pseudofermions and related vectors
  R = malloc(sizeof *R * info.N);
  PF = malloc(sizeof *PF * info.N);
  chi = malloc(sizeof *chi * info.N);

  // Set up neighbors list 'hop'
  // Initialize inverse potential V^{-1} and fermion operator D
  initHop(info, hop);
  initD(info, D, hop);
  initV(info, VInv, hop);
  calcVInv(info, VInv);         // Invert V using LAPACK

  // Initialize random number generator and phi field
  rlxd_init(1, seed);
  ranlxd(phi, info.N);          // Random initialization

  // Generate complex gaussian vector R and use it to
  // set up pseudofermion PF = D.R (adding phi-dependent contributions)
  N = info.N;
  gbt = info.g * info.beta / info.t;    // Relevant combination
  initR(N, R);
  for (i = 0; i < N; i++) {
    PF[i] = 0;
    for (j = 0; j < N; j++)
      PF[i] += D[i][j] * R[j];

    j = (i + N - m) % N;
    PF[i] -= I * gbt * phi[j] * R[j] * (2 * floor((i + N - m) / N) - 1);
  }

  // Reset timer for CG
  dtime += (double)clock();
  dtime /= CLOCKS_PER_SEC;
  printf("\nInitialization done in %.2g sec\n", dtime);
  dtime = -1.0 * (double)clock();

  // Compute chi = (D.Ddag)^{-1} Phi via CG
  printf("\nStarting CG inversion\n");
  CG(info, phi, D, PF, chi, hop);

  dtime += (double)clock();
  dtime /= CLOCKS_PER_SEC;
  printf("\nCG inversion done in %.2g sec\n", dtime);

  return 0;
}
// -----------------------------------------------------------------
