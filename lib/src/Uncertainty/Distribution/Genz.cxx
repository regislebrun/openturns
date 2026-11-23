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
#include <cmath>
#include <algorithm>

#include "openturns/Genz.hxx"
#include "openturns/Exception.hxx"
#include "openturns/SpecFunc.hxx"
#include "openturns/DistFunc.hxx"
#include "openturns/LowDiscrepancySequence.hxx"
#include "openturns/SobolSequence.hxx"
#include "openturns/HaltonSequence.hxx"
#include "openturns/ReverseHaltonSequence.hxx"
#include "openturns/FaureSequence.hxx"
#include "openturns/HaselgroveSequence.hxx"
#include "openturns/TBBImplementation.hxx"
#include "openturns/TriangularMatrix.hxx"

BEGIN_NAMESPACE_OPENTURNS

namespace {

/* Neumaier compensated summation: the weighted estimates add millions of
   terms of very different magnitudes, and the naive sum loses the small ones */
struct CompensatedSum
{
  Scalar total = 0.0;
  Scalar correction = 0.0;

  void add(const Scalar value)
  {
    const Scalar newTotal = total + value;
    if (std::abs(total) >= std::abs(value))
      correction += (total - newTotal) + value;
    else
      correction += (value - newTotal) + total;
    total = newTotal;
  }

  Scalar getValue() const
  {
    return total + correction;
  }

  void join(const CompensatedSum & other)
  {
    add(other.getValue());
  }
};

/* Map a probability to a finite standard normal quantile, keeping it strictly
   inside (0, 1): without the clipping a uniform point equal to 0 on a
   semi-infinite conditional interval would produce an infinite quantile */
Scalar getGenzQuantile(const Scalar probability)
{
  return DistFunc::qNormal(std::min(1.0 - SpecFunc::Precision,
                                    std::max(SpecFunc::Precision, probability)));
}

/* Check the arguments shared by the rectangular probability, moment and
   sampling routines */
void checkGenzArguments(const Point & a,
                        const Point & b,
                        const Point & mu,
                        const TriangularMatrix & L)
{
  const UnsignedInteger d = L.getDimension();
  if (d == 0) throw InvalidDimensionException(HERE) << "Error: the dimension must be at least 1, here d=" << d;
  if (a.getDimension() != d) throw InvalidDimensionException(HERE) << "Error: the lower bound a must have dimension=" << d << ", here dimension=" << a.getDimension();
  if (b.getDimension() != d) throw InvalidDimensionException(HERE) << "Error: the upper bound b must have dimension=" << d << ", here dimension=" << b.getDimension();
  if (mu.getDimension() != d) throw InvalidDimensionException(HERE) << "Error: the mean mu must have dimension=" << d << ", here dimension=" << mu.getDimension();
  for (UnsignedInteger i = 0; i < d; ++i)
  {
    if (!(a[i] <= b[i])) throw InvalidArgumentException(HERE) << "Error: the lower bound a must be less or equal to the upper bound b, here a[" << i << "]=" << a[i] << " > b[" << i << "]=" << b[i];
  }
  for (UnsignedInteger i = 0; i < d; ++i)
  {
    if (!(L(i, i) > 0.0)) throw InvalidArgumentException(HERE) << "Error: the diagonal of the Cholesky factor must be strictly positive, here L(" << i << "," << i << ")=" << L(i, i);
  }
}

/* Replace the infinite bounds by finite values of the same effect. The sign is
   kept: a lower bound equal to +inf, or an upper bound equal to -inf, defines
   an empty interval and must not be turned into a whole line. */
void getFiniteGenzBounds(const Point & a,
                         const Point & b,
                         Point & lower,
                         Point & upper)
{
  lower = a;
  upper = b;
  const Scalar maxScalar = SpecFunc::MaxScalar;
  for (UnsignedInteger i = 0; i < a.getDimension(); ++i)
  {
    if (std::isinf(lower[i])) lower[i] = (lower[i] > 0.0) ? maxScalar : -maxScalar;
    if (std::isinf(upper[i])) upper[i] = (upper[i] > 0.0) ? maxScalar : -maxScalar;
  }
}

/* Center of the second-moment accumulation: the projection of the mean onto
   the box. Centering there instead of at the mean keeps the sampled
   deviations of the order of the box size when the box lies far from the
   mean, which avoids the cancellation in E[X X^T] - E[X] E[X]^T */
Point getMomentCenter(const Point & mu,
                      const Point & lower,
                      const Point & upper)
{
  const UnsignedInteger d = mu.getDimension();
  Point center(d);
  for (UnsignedInteger i = 0; i < d; ++i)
    center[i] = std::min(upper[i], std::max(lower[i], mu[i]));
  return center;
}

/* Sequential Genz decomposition shared by the probability and moment
   entry points. Given the uniform point u, the coordinates of
   z = L^{-1} (x - mu) are drawn in turn from their conditional truncated
   standard normal distribution and p = prod_i (f_i - e_i) is the conditional
   probability of the box, so that P(box) = E[p] and E[g(x) 1_box] = E[p g(x)].
   When WithPoint is true the difference x - mu between the sampled point and
   the mean is stored in meanDeviation and the last coordinate is drawn as
   well. */
template <Bool WithPoint>
Scalar mvn_genz(
    const Point & a,
    const Point & b,
    const TriangularMatrix & L,
    const Sample & sobol,
    const UnsignedInteger k,
    const UnsignedInteger offset,
    const UnsignedInteger d,
    Point & meanDeviation)
{
  Scalar p = 1.0;
  Point y(d, 0.0);
  for (UnsignedInteger i = 0; i < d; ++i)
  {
    // Degenerate dimension: the calling code should have marginalized it
    if (!(L(i, i) > 0.0)) return 0.0;
    const Scalar ai = (a[i] - y[i]) / L(i, i);
    const Scalar bi = (b[i] - y[i]) / L(i, i);
    const Scalar ei = DistFunc::pNormal(ai);
    const Scalar fi = DistFunc::pNormal(bi);
    p *= (fi - ei);
    if (p <= 0.0) break;
    if (WithPoint || (i < d - 1))
    {
      const Scalar zi = getGenzQuantile(ei + sobol(k, i + offset) * (fi - ei));
      if constexpr (WithPoint)
      {
        meanDeviation[i] = y[i] + L(i, i) * zi;
      }
      for (UnsignedInteger j = i + 1; j < d; ++j)
      {
        y[j] += L(j, i) * zi;
      }
    }
  }
  return SpecFunc::Clip01(p);
}

} // anonymous namespace

namespace
{
LowDiscrepancySequence GetGenzSequence(const UnsignedInteger dimension)
{
  const String name(ResourceMap::GetAsString("Genz-LowDiscrepancySequence"));
  if (name == "HaltonSequence") return HaltonSequence(dimension);
  if (name == "ReverseHaltonSequence") return ReverseHaltonSequence(dimension);
  if (name == "FaureSequence") return FaureSequence(dimension);
  if (name == "HaselgroveSequence") return HaselgroveSequence(dimension);
  // default: SobolSequence, also handles "SobolSequence", with HaltonSequence as a fallback if the dimension is too large for SobolSequence
  if (dimension <= SobolSequence::MaximumDimension)
    return SobolSequence(dimension, "MULTIDIGIT");
  return HaltonSequence(dimension, "RANDOM");
}
} // anonymous namespace

Scalar mvn_rectangular_probability(
    const Point & a,
    const Point & b,
    const Point & mu,
    const TriangularMatrix & L,
    const UnsignedInteger n)
{
  const UnsignedInteger d = L.getDimension();

  checkGenzArguments(a, b, mu, L);
  if (n < 1) throw InvalidArgumentException(HERE) << "Error: the sample size n must be positive, here n=" << n;

  const Point a0(a - mu);
  const Point b0(b - mu);

  if (d == 1)
  {
    return DistFunc::pNormal(b0[0] / L(0, 0)) -
           DistFunc::pNormal(a0[0] / L(0, 0));
  }

  const LowDiscrepancySequence seq(GetGenzSequence(d - 1));
  const Sample sobolPoints(seq.generate(n));

  struct Reduce
  {
    const Sample & sobolPoints_;
    const Point & a0_;
    const Point & b0_;
    const TriangularMatrix & L_;
    const UnsignedInteger d_;
    Point unusedPoint_;
    Scalar sum_;

    Reduce(const Sample & sp,
           const Point & a0,
           const Point & b0,
           const TriangularMatrix & L,
           const UnsignedInteger d)
      : sobolPoints_(sp)
      , a0_(a0)
      , b0_(b0)
      , L_(L)
      , d_(d)
      , sum_(0.0)
    {
      // Nothing to do
    }

    Reduce(const Reduce & other,
           TBBImplementation::Split)
      : sobolPoints_(other.sobolPoints_)
      , a0_(other.a0_)
      , b0_(other.b0_)
      , L_(other.L_)
      , d_(other.d_)
      , sum_(0.0)
    {
      // Nothing to do
    }

    void operator()(const TBBImplementation::BlockedRange<UnsignedInteger> & r)
    {
      for (UnsignedInteger k = r.begin(); k != r.end(); ++k)
      {
        sum_ += mvn_genz<false>(a0_, b0_, L_, sobolPoints_, k, 0, d_, unusedPoint_);
      }
    }

    void join(const Reduce& other) { sum_ += other.sum_; }
  };

  Reduce body(sobolPoints, a0, b0, L, d);
  TBBImplementation::ParallelReduce(0, n, body);
  return SpecFunc::Clip01(body.sum_ / n);
}

Scalar mvn_rectangular_probability(
    const Point & a,
    const Point & b,
    const TriangularMatrix & L,
    const UnsignedInteger n)
{
  const UnsignedInteger d = L.getDimension();
  return mvn_rectangular_probability(a, b, Point(d, 0.0), L, n);
}

GenzMoments mvn_rectangular_moments(
    const Point & a,
    const Point & b,
    const Point & mu,
    const TriangularMatrix & L,
    const UnsignedInteger n)
{
  const UnsignedInteger d = L.getDimension();

  checkGenzArguments(a, b, mu, L);
  if (n < 1) throw InvalidArgumentException(HERE) << "Error: the sample size n must be positive, here n=" << n;

  Point lower;
  Point upper;
  getFiniteGenzBounds(a, b, lower, upper);
  const Point a0(lower - mu);
  const Point b0(upper - mu);

  GenzMoments moments;
  moments.probability = 1.0;
  moments.mean = mu;
  moments.covariance = CovarianceMatrix(d);
  moments.effectiveSampleSize = 1.0 * n;

  // Exact one-dimensional formulas: phi and Phi denote the standard normal
  // pdf and cdf, and for a < Z < b
  // E[Z] = (phi(a) - phi(b)) / (Phi(b) - Phi(a)),
  // E[Z^2] = 1 + (a phi(a) - b phi(b)) / (Phi(b) - Phi(a)).
  if (d == 1)
  {
    const Scalar sigma = L(0, 0);
    const Scalar alpha = a0[0] / sigma;
    const Scalar beta = b0[0] / sigma;
    const Scalar probability = DistFunc::pNormal(beta) - DistFunc::pNormal(alpha);
    if (!(probability > 0.0))
      throw InvalidArgumentException(HERE) << "Error: the rectangular probability is zero, the moments are undefined";
    const Scalar alphaPDF = std::isfinite(alpha) ? alpha * DistFunc::dNormal(alpha) : 0.0;
    const Scalar betaPDF = std::isfinite(beta) ? beta * DistFunc::dNormal(beta) : 0.0;
    const Scalar standardMean = (DistFunc::dNormal(alpha) - DistFunc::dNormal(beta)) / probability;
    const Scalar standardSecondMoment = 1.0 + (alphaPDF - betaPDF) / probability;
    moments.probability = SpecFunc::Clip01(probability);
    moments.mean = Point(1, mu[0] + sigma * standardMean);
    moments.covariance(0, 0) = sigma * sigma * std::max(0.0, standardSecondMoment - standardMean * standardMean);
    return moments;
  }

  const LowDiscrepancySequence seq(GetGenzSequence(d));
  const Sample sobolPoints(seq.generate(n));

  const Point center(getMomentCenter(mu, lower, upper));
  const Point shift(mu - center);

  struct MomentsReduce
  {
    const Sample & sobolPoints_;
    const Point & a0_;
    const Point & b0_;
    const Point & shift_;
    const TriangularMatrix & L_;
    const UnsignedInteger d_;
    CompensatedSum probability_;
    CompensatedSum squaredProbability_;
    Collection<CompensatedSum> firstMoment_;
    Collection<CompensatedSum> secondMoment_;
    Point meanDeviation_;
    Point centerDeviation_;

    MomentsReduce(const Sample & sp,
                  const Point & a0,
                  const Point & b0,
                  const Point & shift,
                  const TriangularMatrix & L,
                  const UnsignedInteger d)
      : sobolPoints_(sp)
      , a0_(a0)
      , b0_(b0)
      , shift_(shift)
      , L_(L)
      , d_(d)
      , firstMoment_(d)
      , secondMoment_(d * (d + 1) / 2)
      , meanDeviation_(d, 0.0)
      , centerDeviation_(d, 0.0)
    {
      // Nothing to do
    }

    MomentsReduce(const MomentsReduce & other,
                  TBBImplementation::Split)
      : sobolPoints_(other.sobolPoints_)
      , a0_(other.a0_)
      , b0_(other.b0_)
      , shift_(other.shift_)
      , L_(other.L_)
      , d_(other.d_)
      , firstMoment_(other.d_)
      , secondMoment_(other.d_ * (other.d_ + 1) / 2)
      , meanDeviation_(other.d_, 0.0)
      , centerDeviation_(other.d_, 0.0)
    {
      // Nothing to do
    }

    void operator()(const TBBImplementation::BlockedRange<UnsignedInteger> & r)
    {
      for (UnsignedInteger k = r.begin(); k != r.end(); ++k)
      {
        const Scalar p = mvn_genz<true>(a0_, b0_, L_, sobolPoints_, k, 0, d_, meanDeviation_);
        // A point with a zero weight may carry an infinite coordinate, and
        // 0 * inf would poison the sums
        if (!(p > 0.0)) continue;
        probability_.add(p);
        squaredProbability_.add(p * p);
        for (UnsignedInteger i = 0; i < d_; ++i)
          centerDeviation_[i] = shift_[i] + meanDeviation_[i];
        for (UnsignedInteger i = 0; i < d_; ++i)
        {
          firstMoment_[i].add(p * centerDeviation_[i]);
          for (UnsignedInteger j = 0; j <= i; ++j)
            secondMoment_[i * (i + 1) / 2 + j].add(p * centerDeviation_[i] * centerDeviation_[j]);
        }
      }
    }

    void join(const MomentsReduce & other)
    {
      probability_.join(other.probability_);
      squaredProbability_.join(other.squaredProbability_);
      for (UnsignedInteger i = 0; i < d_; ++i)
        firstMoment_[i].join(other.firstMoment_[i]);
      for (UnsignedInteger i = 0; i < secondMoment_.getSize(); ++i)
        secondMoment_[i].join(other.secondMoment_[i]);
    }
  };

  MomentsReduce body(sobolPoints, a0, b0, shift, L, d);
  TBBImplementation::ParallelReduce(0, n, body);

  const Scalar sumProbability = body.probability_.getValue();
  const Scalar sumSquaredProbability = body.squaredProbability_.getValue();
  // The truncated moments are undefined when the box has a zero probability;
  // this also catches an underflowed weight in a very small box
  if (!(sumProbability > 0.0))
    throw InvalidArgumentException(HERE) << "Error: the rectangular probability is zero, the moments are undefined; increase Genz-DefaultSampleSize";
  moments.probability = SpecFunc::Clip01(sumProbability / n);
  moments.effectiveSampleSize = sumProbability * sumProbability / sumSquaredProbability;
  Point centeredMean(d);
  for (UnsignedInteger i = 0; i < d; ++i)
    centeredMean[i] = body.firstMoment_[i].getValue() / sumProbability;
  moments.mean = Point(d);
  for (UnsignedInteger i = 0; i < d; ++i)
    moments.mean[i] = center[i] + centeredMean[i];
  SymmetricMatrix covariance(d);
  for (UnsignedInteger i = 0; i < d; ++i)
  {
    for (UnsignedInteger j = 0; j <= i; ++j)
    {
      Scalar value = body.secondMoment_[i * (i + 1) / 2 + j].getValue() / sumProbability
                     - centeredMean[i] * centeredMean[j];
      // A round-off negative variance is meaningless
      if (i == j) value = std::max(0.0, value);
      covariance(i, j) = value;
    }
  }
  moments.covariance = CovarianceMatrix(covariance);
  return moments;
}

Scalar mvt_rectangular_probability(
    const Point & a,
    const Point & b,
    const Point & mu,
    const TriangularMatrix & L,
    const Scalar nu,
    const UnsignedInteger n)
{
  const UnsignedInteger d = L.getDimension();

  checkGenzArguments(a, b, mu, L);
  if (n < 1) throw InvalidArgumentException(HERE) << "Error: the sample size n must be positive, here n=" << n;
  if (!(nu > 0.0)) throw InvalidArgumentException(HERE) << "Error: the degrees of freedom nu must be positive, here nu=" << nu;

  const Point a0(a - mu);
  const Point b0(b - mu);

  const LowDiscrepancySequence seq(GetGenzSequence(d));
  const Sample sobolPoints(seq.generate(n));

  struct Reduce
  {
    const Sample & sobolPoints_;
    const Point & a0_;
    const Point & b0_;
    const TriangularMatrix & L_;
    const UnsignedInteger d_;
    const Scalar nu_;
    const Scalar isqrtnu_;
    Point unusedPoint_;
    Scalar sum_;

    Reduce(const Sample & sp,
           const Point & a0,
           const Point & b0,
           const TriangularMatrix & L,
           const UnsignedInteger d,
           const Scalar nu)
      : sobolPoints_(sp)
      , a0_(a0)
      , b0_(b0)
      , L_(L)
      , d_(d)
      , nu_(nu)
      , isqrtnu_(1.0 / std::sqrt(nu))
      , sum_(0.0)
    {
      // Nothing to do
    }

    Reduce(const Reduce & other,
           TBBImplementation::Split)
      : sobolPoints_(other.sobolPoints_)
      , a0_(other.a0_)
      , b0_(other.b0_)
      , L_(other.L_)
      , d_(other.d_)
      , nu_(other.nu_)
      , isqrtnu_(other.isqrtnu_)
      , sum_(0.0)
    {
      // Nothing to do
    }

    void operator()(const TBBImplementation::BlockedRange<UnsignedInteger> & r)
    {
      for (UnsignedInteger k = r.begin(); k != r.end(); ++k)
      {
        const Scalar s = DistFunc::qChi(nu_, sobolPoints_(k, 0)) * isqrtnu_;
        sum_ += mvn_genz<false>(a0_ * s, b0_ * s, L_, sobolPoints_, k, 1, d_, unusedPoint_);
      }
    }

    void join(const Reduce& other) { sum_ += other.sum_; }
  };

  Reduce body(sobolPoints, a0, b0, L, d, nu);
  TBBImplementation::ParallelReduce(0, n, body);
  return SpecFunc::Clip01(body.sum_ / n);
}

Scalar mvt_rectangular_probability(
    const Point & a,
    const Point & b,
    const TriangularMatrix & L,
    const Scalar nu,
    const UnsignedInteger n)
{
  const UnsignedInteger d = L.getDimension();
  return mvt_rectangular_probability(a, b, Point(d, 0.0), L, nu, n);
}

END_NAMESPACE_OPENTURNS
