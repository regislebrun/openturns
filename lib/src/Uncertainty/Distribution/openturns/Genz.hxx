//                                               -*- C++ -*-
/**
 *  @brief Genz algorithm to compute rectangular probabilities of multivariate
 *         Normal and Student distributions, and moments of truncated
 *         multivariate Normal distributions.
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
#ifndef OPENTURNS_GENZ_HXX
#define OPENTURNS_GENZ_HXX

#include "openturns/OTprivate.hxx"
#include "openturns/ResourceMap.hxx"
#include "openturns/Point.hxx"
#include "openturns/CovarianceMatrix.hxx"
#include "openturns/TriangularMatrix.hxx"
BEGIN_NAMESPACE_OPENTURNS

/* Moments of the normal distribution N(mu, L L^T) truncated to the box
   [a, b]. They are exact in dimension 1, and estimated otherwise by the
   weighted Genz quasi-Monte Carlo estimator, which draws the truncated point
   and accumulates its box weight in a single pass. Degenerate coordinates
   (zero diagonal entry of L) must be removed by the caller, as done by
   Normal::computeProbability. */
struct GenzMoments
{
  /* Probability of the box */
  Scalar probability = 0.0;

  /* Mean of the truncated distribution */
  Point mean;

  /* Covariance of the truncated distribution */
  CovarianceMatrix covariance;

  /* Effective sample size (sum p)^2 / sum p^2 of the weighted estimator, a
     tail-degeneracy diagnostic: the estimate is unreliable when it is small
     compared with the number of quasi-Monte Carlo points. It equals the
     number of points when the moments are computed exactly, in dimension 1 */
  Scalar effectiveSampleSize = 0.0;
};

OT_API Scalar mvn_rectangular_probability(
    const Point & a,
    const Point & b,
    const Point & mu,
    const TriangularMatrix & L,
    const UnsignedInteger n = ResourceMap::GetAsUnsignedInteger("Genz-DefaultSampleSize"));

OT_API Scalar mvn_rectangular_probability(
    const Point & a,
    const Point & b,
    const TriangularMatrix & L,
    const UnsignedInteger n = ResourceMap::GetAsUnsignedInteger("Genz-DefaultSampleSize"));

OT_API Scalar mvt_rectangular_probability(
    const Point & a,
    const Point & b,
    const Point & mu,
    const TriangularMatrix & L,
    const Scalar nu,
    const UnsignedInteger n = ResourceMap::GetAsUnsignedInteger("Genz-DefaultSampleSize"));

OT_API Scalar mvt_rectangular_probability(
    const Point & a,
    const Point & b,
    const TriangularMatrix & L,
    const Scalar nu,
    const UnsignedInteger n = ResourceMap::GetAsUnsignedInteger("Genz-DefaultSampleSize"));

/* Moments of the normal distribution N(mu, L L^T) truncated to [a, b]. The
   infinite bounds are handled as in mvn_rectangular_probability, and the
   moment and squared probability sums use compensated summation. A box with
   a zero probability is reported as an InvalidArgumentException. */
OT_API GenzMoments mvn_rectangular_moments(
    const Point & a,
    const Point & b,
    const Point & mu,
    const TriangularMatrix & L,
    const UnsignedInteger n = ResourceMap::GetAsUnsignedInteger("Genz-DefaultSampleSize"));

END_NAMESPACE_OPENTURNS

#endif // OPENTURNS_GENZ_HXX
