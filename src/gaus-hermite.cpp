/*
  Copyright (c) 2014, Alexander W Blocker

  Permission is hereby granted, free of charge, to any person obtaining
  a copy of this software and associated documentation files (the
  "Software"), to deal in the Software without restriction, including
  without limitation the rights to use, copy, modify, merge, publish,
  distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to
  the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/* The code is modified by Benjamin Christoffersen */

#include "gaus-hermite.h"
#include <cmath>
#include <math.h>

#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCLEN
#define FCLEN
#endif
#ifndef FCONE
#define FCONE
#endif
#include <R_ext/Lapack.h>

namespace GaussHermite {
using std::vector;
using std::abs;

void buildHermiteJacobi(vector<double> &D, vector<double> &E){
  // Construct symmetric tridiagonal matrix similar to Jacobi matrix
  // for Hermite polynomials
  //
  // On exit, D contains diagonal elements of said matrix;
  // E contains subdiagonal elements.
  //
  // Need D of size n, E of size n-1
  //
  // Building matrix based on recursion relation for monic versions of Hermite
  // polynomials:
  //      p_n(x) = H_n(x) / 2^n
  //      p_n+1(x) + (B_n-x)*p_n(x) + A_n*p_n-1(x) = 0
  //      B_n = 0
  //      A_n = n/2
  //
  // Matrix similar to Jacobi (J) defined by:
  //      J_i,i = B_i-1, i = 1, ..., n
  //      J_i,i-1 = J_i-1,i = sqrt(A_i-1), i = 2, ..., n
  unsigned const n = D.size();
  for(unsigned i = 0; i < n; i++){
    // Build diagonal
    D[i] = 0.;
    // Build sub/super-diagonal
    E[i] = std::sqrt((i + 1.) / 2.);
  }
}

void quadInfoGolubWelsch
  (vector<double> &D, vector<double> &E, double const mu0,
   vector<double> &x, vector<double> &w) {
  // Compute weights & nodes for Gaussian quadrature using Golub-Welsch
  // algorithm.
  //
  // First need to build symmetric tridiagonal matrix J similar to Jacobi for
  // desired orthogonal polynomial (based on recurrence relation).
  //
  // D contains the diagonal of this matrix J, and E contains the
  // sub/super-diagonal.
  //
  // This routine finds the eigenvectors & values of the given J matrix.
  //
  // The eigenvalues correspond to the nodes for the desired quadrature rule
  // (roots of the orthogonal polynomial).
  //
  // The eigenvectors can be used to compute the weights for the quadrature rule
  // via:
  //
  //      w_j = mu0 * (v_j,1)^2
  //
  // where mu0 = \int_a^b w(x) dx
  // (integral over range of integration of weight function)
  //
  // and
  //
  // v_j,1 is the first entry of the jth normalized (to unit length)
  // eigenvector.
  //
  // On exit, x (length n) contains nodes for quadrature rule, and w (length n)
  // contains weights for quadrature rule.
  //
  // Note that contents of D & E are destroyed on exit
  //

  // Setup for eigenvalue computations
  char JOBZ = 'V';  // Flag to compute both eigenvalues & vectors.
  int INFO;
  int const n = x.size();
  vector<double> WORK(2 * n - 2);
  vector<double> Z(n * n);  // This holds the resulting eigenvectors.

  // Run eigen decomposition
  F77_NAME(dstev)(&JOBZ, &n, &D[0], &E[0],  // Job flag & input matrix
                  &Z[0], &n,       // Output array for eigenvectors & dim
                  &WORK[0], &INFO  // Workspace & info flag
                  FCONE);

  // Setup x & w
  for (int i = 0; i < n; i++) {
    x[i] = D[i];
    w[i] = mu0 * Z[i * n] * Z[i * n];
  }
}

int gaussHermiteDataGolubWelsch(vector<double> &x, vector<double> &w) {
  //
  // Calculates nodes & weights for Gauss-Hermite integration of order n
  //
  // Need x & w of same size n
  //
  // Using standard formulation (no generalizations or polynomial adjustment)
  //
  // Evaluations use Golub-Welsch algorithm; numerically stable for n>=100
  //
  // Build Jacobi-similar symmetric tridiagonal matrix via diagonal &
  // sub-diagonal
  unsigned const n = x.size();
  vector<double> D(n), E(n);
  buildHermiteJacobi(D, E);

  // Get nodes & weights
  double const mu0 = sqrt(M_PI);
  quadInfoGolubWelsch(D, E, mu0, x, w);

  return 0;
}

HermiteData<double> GaussHermiteData(unsigned const n) {
  HermiteData<double> out(n);

  // Build Gauss-Hermite integration rules
  gaussHermiteDataGolubWelsch(out.x, out.w);

  return out;
}

} // namespace GaussHermite
