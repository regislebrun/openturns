//                                               -*- C++ -*-
/**
 *  @brief The Rice distribution
 *
 *  Copyright 2005-2025 Airbus-EDF-IMACS-ONERA-Phimeca
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
#include "openturns/Rice.hxx"
#include "openturns/SpecFunc.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/DistFunc.hxx"
#include "openturns/SpecFunc.hxx"
#include "openturns/Distribution.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(Rice)

static const Factory<Rice> Factory_Rice;

/* Default constructor */
Rice::Rice()
  : DistributionImplementation()
  , beta_(-1.0)
  , nu_(-1.0)
  , maximumIteration_(ResourceMap::GetAsUnsignedInteger("DistFunc-MaximumIteration"))
{
  setName("Rice");
  setNu(0.0);
  setBeta(1.0);
  setDimension(1);
  computeRange();
}

/* Parameters constructor */
Rice::Rice(const Scalar beta,
           const Scalar nu)
  : DistributionImplementation()
  , beta_(0.0)
  , nu_(nu)
  , maximumIteration_(ResourceMap::GetAsUnsignedInteger("DistFunc-MaximumIteration"))
{
  setName("Rice");
  setNu(nu);
  // This call sets also the range
  setBeta(beta);
  setDimension(1);
}

/* Comparison operator */
Bool Rice::operator ==(const Rice & other) const
{
  if (this == &other) return true;
  return (beta_ == other.beta_) && (nu_ == other.nu_);
}

Bool Rice::equals(const DistributionImplementation & other) const
{
  const Rice* p_other = dynamic_cast<const Rice*>(&other);
  return p_other && (*this == *p_other);
}

/* String converter */
String Rice::__repr__() const
{
  OSS oss;
  oss << "class=" << Rice::GetClassName()
      << " name=" << getName()
      << " dimension=" << getDimension()
      << " beta=" << beta_
      << " nu=" << nu_;
  return oss;
}

String Rice::__str__(const String & ) const
{
  OSS oss;
  oss << getClassName() << "(beta = " << beta_ << ", nu = " << nu_ << ")";
  return oss;
}

/* Virtual constructor */
Rice * Rice::clone() const
{
  return new Rice(*this);
}

/* Compute the numerical range of the distribution given the parameters values */
void Rice::computeRange()
{
  const Point lowerBound(1, 0.0 );
  const Point upperBound(1, computeScalarQuantile(cdfEpsilon_, true));
  const Interval::BoolCollection finiteLowerBound(1, true);
  const Interval::BoolCollection finiteUpperBound(1, false);
  setRange(Interval(lowerBound, upperBound, finiteLowerBound, finiteUpperBound));
}


/* Get one realization of the distribution */
Point Rice::getRealization() const
{
  const Scalar x = beta_ * DistFunc::rNormal() + nu_;
  const Scalar y = beta_ * DistFunc::rNormal();
  return Point(1.0, std::sqrt(x * x + y * y));
}


/* Get the PDF of the distribution */
Scalar Rice::computePDF(const Point & point) const
{
  if (point.getDimension() != 1) throw InvalidArgumentException(HERE) << "Error: the given point must have dimension=1, here dimension=" << point.getDimension();

  const Scalar x = point[0];
  if (x <= 0.0) return 0.0;
  const Scalar xScaled = x / beta_;
  const Scalar nuScaled = nu_ / beta_;
  return xScaled / beta_ * std::exp(-0.5 * (xScaled * xScaled + nuScaled * nuScaled) + SpecFunc::LogBesselI0(xScaled * nuScaled));
}


/* Get the logarithm of the PDF of the distribution */
Scalar Rice::computeLogPDF(const Point & point) const
{
  if (point.getDimension() != 1) throw InvalidArgumentException(HERE) << "Error: the given point must have dimension=1, here dimension=" << point.getDimension();

  const Scalar x = point[0];
  if (x <= 0.0) return SpecFunc::LowestScalar;
  const Scalar xScaled = x / beta_;
  const Scalar nuScaled = nu_ / beta_;
  return std::log(xScaled / beta_) - 0.5 * (xScaled * xScaled + nuScaled * nuScaled) + SpecFunc::LogBesselI0(xScaled * nuScaled);
}


/* Get the CDF of the distribution */
Scalar Rice::computeCDF(const Point & point) const
{
  if (point.getDimension() != 1) throw InvalidArgumentException(HERE) << "Error: the given point must have dimension=1, here dimension=" << point.getDimension();

  if (point[0] <= 0.0) return 0.0;
  const Scalar lambda = std::pow(nu_ / beta_, 2);
  const Scalar y = std::pow(point[0] / beta_, 2);
  return DistFunc::pNonCentralChiSquare(2, lambda, y, false, pdfEpsilon_, maximumIteration_);
}

Scalar Rice::computeComplementaryCDF(const Point & point) const
{
  if (point.getDimension() != 1) throw InvalidArgumentException(HERE) << "Error: the given point must have dimension=1, here dimension=" << point.getDimension();

  if (point[0] <= 0.0) return 1.0;
  const Scalar lambda = std::pow(nu_ / beta_, 2);
  const Scalar y = std::pow(point[0] / beta_, 2);
  return DistFunc::pNonCentralChiSquare(2, lambda, y, true, pdfEpsilon_, maximumIteration_);
}

Scalar Rice::computeProbability(const Interval & interval) const
{
  if (interval.getDimension() != 1)
    throw InvalidArgumentException(HERE) << "computeProbability expected an interval of dimension=" << dimension_ << ", got dimension=" << interval.getDimension();
  return computeProbabilityGeneral1D(interval.getLowerBound()[0], interval.getUpperBound()[0]);
}

/* Compute the mean of the distribution */
void Rice::computeMean() const
{
  // 1.253314137315500251207882 = sqrt(pi/2)
  //const Scalar x = -0.5 * std::pow(nu_ / beta_, 2);
  //mean_ = Point(1, beta_ * 1.253314137315500251207882 * SpecFunc::HyperGeom_1_1(-0.5, 1, x));
  const Scalar t = std::pow(0.5 * nu_ / beta_, 2);
  const Scalar B0 = SpecFunc::BesselI0(t);
  const Scalar B1 = SpecFunc::BesselI1(t);
  const Scalar E = std::exp(-t);
  const Scalar L = ((1.0 + 2.0 * t) * B0 + 2.0 * t * B1) * E;
  mean_ = Point(1, beta_ * 1.253314137315500251207882 * L);
}

/* Get the standard deviation of the distribution */
Point Rice::getStandardDeviation() const
{
  if (!isAlreadyComputedCovariance_) computeCovariance();
  return Point(1, std::sqrt(covariance_(0, 0)));
}

/* Get the skewness of the distribution */
Point Rice::getSkewness() const
{
  // 1.253314137315500251207882 = sqrt(pi/2)
  const Scalar t = std::pow(0.5 * nu_ / beta_, 2);
  const Scalar B0 = SpecFunc::BesselI0(t);
  const Scalar B1 = SpecFunc::BesselI1(t);
  const Scalar E = std::exp(-t);
  const Scalar L = ((1.0 + 2.0 * t) * B0 + 2.0 * t * B1) * E;
  const Scalar mu3 = 1.253314137315500251207882 * beta_ * (beta_ * beta_ * (M_PI * std::pow(L, 3) - B0 * E - 2.0 * L) - 2.0 * L * std::pow(nu_, 2));
  const Scalar sigma3 = std::pow(getStandardDeviation()[0], 3);
  return Point(1, mu3 / sigma3);
}

/* Get the kurtosis of the distribution */
Point Rice::getKurtosis() const
{
  // 1.253314137315500251207882 = sqrt(pi/2)
  const Scalar nu2 = nu_ * nu_;
  const Scalar beta2 = beta_ * beta_;
  const Scalar t = 0.25 * nu2 / beta2;
  const Scalar B0 = SpecFunc::BesselI0(t);
  const Scalar B1 = SpecFunc::BesselI1(t);
  const Scalar E = std::exp(-t);
  const Scalar L = ((1.0 + 2.0 * t) * B0 + 2.0 * t * B1) * E;
  const Scalar L2 = L * L;
  const Scalar mu4 = nu2 * nu2 + beta2 * (nu2 * (8.0 + M_PI * L2) + beta2 * (8.0 + 2 * M_PI * E * L * B0 - 2.0 * M_PI * L2 - 0.75 * M_PI * M_PI * L2 * L2));
  const Scalar sigma4 = std::pow(getStandardDeviation()[0], 4);
  return Point(1, mu4 / sigma4);
}

/* Get the standard representative in the parametric family, associated with the standard moments */
Distribution Rice::getStandardRepresentative() const
{
  Rice standard(1.0, 0.0);
  standard.setDescription(getDescription());
  return standard;
}


/* Compute the covariance of the distribution */
void Rice::computeCovariance() const
{
  covariance_ = CovarianceMatrix(1);
  Scalar covariance = 0.0;
  const Scalar mu = getMean()[0];
  covariance = 2.0 * beta_ * beta_ + (nu_ - mu) * (nu_ + mu);
  covariance_(0, 0) = covariance;
  isAlreadyComputedCovariance_ = true;
}

/* Parameters value accessor */
Point Rice::getParameter() const
{
  return {beta_, nu_};
}

void Rice::setParameter(const Point & parameter)
{
  if (parameter.getSize() != 2) throw InvalidArgumentException(HERE) << "Error: expected 2 values, got " << parameter.getSize();
  const Scalar w = getWeight();
  *this = Rice(parameter[0], parameter[1]);
  setWeight(w);
}

/* Parameters description accessor */
Description Rice::getParameterDescription() const
{
  return {"beta", "nu"};
}

/* Beta accessor */
void Rice::setBeta(const Scalar beta)
{
  if (!(beta > 0.0)) throw InvalidArgumentException(HERE) << "Beta MUST be positive";
  if (beta != beta_)
  {
    beta_ = beta;
    isAlreadyComputedMean_ = false;
    isAlreadyComputedCovariance_ = false;
    computeRange();
  }
}

Scalar Rice::getBeta() const
{
  return beta_;
}

/* Nu accessor */
void Rice::setNu(const Scalar nu)
{
  if (!(nu >= 0.0)) throw InvalidArgumentException(HERE) << "Nu MUST be positive";
  if (nu != nu_)
  {
    nu_ = nu;
    isAlreadyComputedMean_ = false;
    // The covariance does not depend on nu
    computeRange();
  }
}

Scalar Rice::getNu() const
{
  return nu_;
}

/* Maximum iterations accessor */
void Rice::setMaximumIteration(const UnsignedInteger maximumIteration)
{
  maximumIteration_ = maximumIteration;
}

UnsignedInteger Rice::getMaximumIteration() const
{
  return maximumIteration_;
}

/* Method save() stores the object through the StorageManager */
void Rice::save(Advocate & adv) const
{
  DistributionImplementation::save(adv);
  adv.saveAttribute( "beta_", beta_ );
  adv.saveAttribute( "nu_", nu_ );
  adv.saveAttribute( "maximumIteration_", maximumIteration_ );
}

/* Method load() reloads the object from the StorageManager */
void Rice::load(Advocate & adv)
{
  DistributionImplementation::load(adv);
  if (adv.hasAttribute("sigma_")) // old parameters
    adv.loadAttribute( "sigma_", beta_ );
  else
    adv.loadAttribute( "beta_", beta_ );
  adv.loadAttribute( "nu_", nu_ );
  adv.loadAttribute( "maximumIteration_", maximumIteration_ );
  computeRange();
}

END_NAMESPACE_OPENTURNS
