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
#include "openturns/FastGaussQuadrature.hxx"
#include "openturns/Point.hxx"
#include "openturns/Lapack.hxx"
#include "openturns/Exception.hxx"
#include <cmath>
#include <algorithm>
#include <limits>

BEGIN_NAMESPACE_OPENTURNS

/**
 * @namespace FastGaussQuadrature
 *
 * Shared polished Gauss quadrature solver.
 *
 * The symmetric three-term recurrence for orthonormal polynomials is:
 *   x * p_j(x) = b_{j+1} * p_{j+1}(x) + gamma_j * p_j(x) + b_j * p_{j-1}(x)
 * with b_0 = 0, p_{-1} = 0, p_0 = 1.
 *
 * The Jacobi matrix J has:
 *   J[i,i]   = gamma[i]       (diagonal)
 *   J[i,i+1] = b[i+1]        (off-diagonal, i = 0..n-2)
 *
 * Nodes are eigenvalues of J; weights are computed from the Christoffel identity
 * or equivalently w_i = 1 / (p_{n-1}(x_i) * p'_n(x_i)).
 */
namespace FastGaussQuadrature
{
  void PolishedSolve(const Scalar * gamma,
                     const Scalar * b,
                     const UnsignedInteger n,
                     Scalar * nodes,
                     Scalar * weights)
  {
    if (n == 0) throw InvalidArgumentException(HERE) << "Error: n must be > 0";
    if (n == 1)
    {
      nodes[0] = gamma[0];
      weights[0] = 1.0;
      return;
    }

    // Copy diagonal and off-diagonal for LAPACK
    Scalar * d = nodes;  // eigenvalues written in-place to nodes
    Point e(n - 1);
    for (UnsignedInteger i = 0; i < n - 1; ++i)
      e[i] = b[i + 1];

    // Eigenvalues only (jobz='N')
    char jobz = 'N';
    int ljobz = 1;
    int lwork = std::max(1, 2 * static_cast<int>(n) - 2);
    Point work(lwork);
    Point zDummy(1);
    int info = 0;
    int size = static_cast<int>(n);

    dstev_(&jobz, &size, &d[0], &e[0], &zDummy[0], &size, &work[0], &info, &ljobz);
    if (info != 0) throw InternalException(HERE) << "LAPACK DSTEV error: info=" << info;

    // Newton refinement of every node
    const Scalar tolerance = 1.0e-15;
    Point fm(n);
    Point dfm(n);
    for (UnsignedInteger i = 0; i < n; ++i)
    {
      Scalar x = d[i];
      Scalar previousMagnitude = 0.0;
      UnsignedInteger stallCount = 0;
      for (UnsignedInteger iter = 0; iter < 20; ++iter)
      {
        // Evaluate p_n(x) and p'_n(x) by the symmetric three-term recurrence
        // Loop runs n-1 times to compute p_{n-1} and p'_{n-1}
        Scalar pp = 0.0;
        Scalar pk = 1.0;
        Scalar dpp = 0.0;
        Scalar dpk = 0.0;
        for (UnsignedInteger j = 0; j < n - 1; ++j)
        {
          const Scalar bjp1 = b[j + 1];
          const Scalar pjp1 = ((x - gamma[j]) * pk - b[j] * pp) / bjp1;
          const Scalar dpjp1 = ((x - gamma[j]) * dpk + pk - b[j] * dpp) / bjp1;
          pp = pk;
          pk = pjp1;
          dpp = dpk;
          dpk = dpjp1;
        }
        // Now pp = p_{n-2}, pk = p_{n-1}, dpk = p'_{n-1}
        // Compute p_n and p'_n using the recurrence at the last step (b[n] = 0)
        const Scalar pN = ((x - gamma[n - 1]) * pk - b[n - 1] * pp);
        const Scalar dPN = ((x - gamma[n - 1]) * dpk + pk - b[n - 1] * dpp);

        const Scalar step = pN / dPN;
        x -= step;
        const Scalar magnitude = std::abs(step) / (1.0 + std::abs(x));
        if (magnitude <= tolerance) break;
        if ((iter > 0) && (magnitude >= 0.5 * previousMagnitude))
        {
          ++stallCount;
          if (stallCount >= 2) break;
        }
        else stallCount = 0;
        previousMagnitude = magnitude;
      }
      // Re-evaluate at the final point for weight computation
      // We need p_{n-1}(x_i) and p'_n(x_i)
      {
        Scalar pp = 0.0;
        Scalar pk = 1.0;
        Scalar dpp = 0.0;
        Scalar dpk = 0.0;
        for (UnsignedInteger j = 0; j < n - 1; ++j)
        {
          const Scalar bjp1 = b[j + 1];
          const Scalar pjp1 = ((x - gamma[j]) * pk - b[j] * pp) / bjp1;
          const Scalar dpjp1 = ((x - gamma[j]) * dpk + pk - b[j] * dpp) / bjp1;
          pp = pk;
          pk = pjp1;
          dpp = dpk;
          dpk = dpjp1;
        }
        // pp = p_{n-2}, pk = p_{n-1}, dpk = p'_{n-1}
        fm[i] = pk;   // p_{n-1}
        dfm[i] = ((x - gamma[n - 1]) * dpk + pk - b[n - 1] * dpp);  // p'_n
      }
      d[i] = x;
    }

    // Weights from 1/(p_{n-1}(x) * p'_n(x)), rescaled to avoid over/underflow
    {
      Scalar low = std::log(std::numeric_limits<Scalar>::max());
      Scalar high = -low;
      for (UnsignedInteger i = 0; i < n; ++i)
      {
        const Scalar magnitude = std::abs(fm[i]);
        if (magnitude == 0.0) continue;
        const Scalar logMagnitude = std::log(magnitude);
        low = std::min(low, logMagnitude);
        high = std::max(high, logMagnitude);
      }
      for (UnsignedInteger i = 0; i < n; ++i)
      {
        const Scalar magnitude = std::abs(dfm[i]);
        if (magnitude == 0.0) continue;
        const Scalar logMagnitude = std::log(magnitude);
        low = std::min(low, logMagnitude);
        high = std::max(high, logMagnitude);
      }
      if (high > low)
      {
        const Scalar shift = std::exp(0.5 * (low + high));
        for (UnsignedInteger i = 0; i < n; ++i)
        {
          fm[i] /= shift;
          dfm[i] /= shift;
        }
      }
    }
    Scalar totalWeight = 0.0;
    for (UnsignedInteger i = 0; i < n; ++i)
    {
      weights[i] = 1.0 / (fm[i] * dfm[i]);
      totalWeight += weights[i];
    }
    // Normalize weights to sum to 1 (pdf convention)
    for (UnsignedInteger i = 0; i < n; ++i)
      weights[i] /= totalWeight;
  }
} // namespace FastGaussQuadrature

END_NAMESPACE_OPENTURNS
