//                                               -*- C++ -*-
/**
 *  @brief Shared polished Gauss quadrature solver
 *
 *  Copyright 2005-2026 Airbus-EDF-IMACS-ONERA-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef OPENTURNS_FASTGAUSSQUADRATURE_HXX
#define OPENTURNS_FASTGAUSSQUADRATURE_HXX

#include "openturns/OTprivate.hxx"

BEGIN_NAMESPACE_OPENTURNS

/**
 * @namespace FastGaussQuadrature
 *
 * Shared polished Gauss quadrature solver used by FastHermite, FastLaguerre
 * and FastJacobi.
 *
 * The symmetric three-term recurrence for orthonormal polynomials is:
 *   x * p_j(x) = b_{j+1} * p_{j+1}(x) + gamma_j * p_j(x) + b_j * p_{j-1}(x)
 * with b_0 = 0, p_{-1} = 0, p_0 = 1.
 *
 * Algorithm (the "polished" method from the report):
 * 1. Build symmetric tridiagonal Jacobi matrix: d[i] = gamma[i], e[i] = b[i+1]
 * 2. LAPACK dstev_ with jobz='N' (eigenvalues only, no eigenvectors)
 * 3. Newton refinement of every node using the three-term recurrence
 * 4. Weights from w_i = 1 / (p_{n-1}(x_i) * p'_n(x_i)), normalized to sum to 1
 *
 * This avoids the eigenvector path entirely, which is the source of accuracy
 * degrades for exponentially decaying weights (Laguerre at large n).
 *
 * Reference: Golub, G. H. and Welsch, J. H. (1969).
 *            Calculation of Gauss Quadrature Rules.
 *            Mathematics of Computation, 23(106), 221-230.
 *
 * Reference: Yakimiw, E. (1996).
 *            Accurate computation of weights of Gauss-Laguerre and
 *            Gauss-Hermite quadrature formulae.
 */
namespace FastGaussQuadrature
{
  /** Polished Gauss quadrature: eigenvalues-only tridiagonal solve + Newton
   *  refinement + Christoffel weight recovery.
   *
   *  @param gamma  diagonal of the Jacobi matrix, length n
   *  @param b      off-diagonal of the Jacobi matrix, length n; b[0] is unused,
   *                b[i] for i=1..n-1 are the sub/super-diagonal entries
   *  @param n      number of quadrature nodes (> 0)
   *  @param nodes  input: initial guesses (length n); output: refined nodes
   *  @param weights output: quadrature weights (length n), sums to 1
   */
  void PolishedSolve(const Scalar * gamma,
                     const Scalar * b,
                     const UnsignedInteger n,
                     Scalar * nodes,
                     Scalar * weights);
} // namespace FastGaussQuadrature

END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_FASTGAUSSQUADRATURE_HXX */
