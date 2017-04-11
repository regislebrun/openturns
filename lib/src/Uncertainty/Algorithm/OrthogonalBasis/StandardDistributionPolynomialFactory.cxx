//                                               -*- C++ -*-
/**
 *  @brief OrthogonalUniVariatePolynomialStandardDistribution polynomial factory
 *
 *  Copyright 2005-2017 Airbus-EDF-IMACS-Phimeca
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
 *  You should have received a copy of the GNU Lesser General Public
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "openturns/StandardDistributionPolynomialFactory.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/Uniform.hxx"
#include "openturns/AdaptiveStieltjesAlgorithm.hxx"
#include "openturns/CharlierFactory.hxx"
#include "openturns/ChebychevFactory.hxx"
#include "openturns/HermiteFactory.hxx"
#include "openturns/JacobiFactory.hxx"
#include "openturns/KrawtchoukFactory.hxx"
#include "openturns/LaguerreFactory.hxx"
#include "openturns/LegendreFactory.hxx"
#include "openturns/HistogramPolynomialFactory.hxx"
#include "openturns/MeixnerFactory.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(StandardDistributionPolynomialFactory);

static const Factory<StandardDistributionPolynomialFactory> Factory_StandardDistributionPolynomialFactory;


/* Default constructor */
StandardDistributionPolynomialFactory::StandardDistributionPolynomialFactory()
  : OrthogonalUniVariatePolynomialFactory(Uniform())
  , orthonormalizationAlgorithm_(AdaptiveStieltjesAlgorithm(Uniform()))
  , specificFamily_()
  , hasSpecificFamily_(false)
  , linear_(1.0)
  , constant_(0.0)
{
  // Initialize the coefficient cache
  initializeCache();
}


/* Parameter constructor */
StandardDistributionPolynomialFactory::StandardDistributionPolynomialFactory(const Distribution & measure)
  : OrthogonalUniVariatePolynomialFactory(measure.getStandardRepresentative())
  , orthonormalizationAlgorithm_(AdaptiveStieltjesAlgorithm(measure.getStandardRepresentative()))
  , specificFamily_()
  , hasSpecificFamily_(false)
  , linear_(1.0)
  , constant_(0.0)
{
  checkSpecificFamily();
  initializeCache();
}


/* Parameter constructor */
StandardDistributionPolynomialFactory::StandardDistributionPolynomialFactory(const OrthonormalizationAlgorithm & orthonormalizationAlgorithm)
  : OrthogonalUniVariatePolynomialFactory(orthonormalizationAlgorithm.getMeasure())
  , orthonormalizationAlgorithm_(orthonormalizationAlgorithm)
  , specificFamily_()
  , hasSpecificFamily_(false)
  , linear_(1.0)
  , constant_(0.0)
{
  checkSpecificFamily();
  initializeCache();
}


/* Virtual constructor */
StandardDistributionPolynomialFactory * StandardDistributionPolynomialFactory::clone() const
{
  return new StandardDistributionPolynomialFactory(*this);
}


/* Calculate the coefficients of recurrence a0n, a1n, a2n such that
   Pn+1(x) = (a0n * x + a1n) * Pn(x) + a2n * Pn-1(x) */
StandardDistributionPolynomialFactory::Coefficients StandardDistributionPolynomialFactory::getRecurrenceCoefficients(const UnsignedInteger n) const
{
  if (hasSpecificFamily_) return specificFamily_.getRecurrenceCoefficients(n);
  else return orthonormalizationAlgorithm_.getRecurrenceCoefficients(n);
}

/* Check the existence of a specific family more efficient for the given measure */
void StandardDistributionPolynomialFactory::checkSpecificFamily()
{
  // Check for special cases. Need a more elegant conception and implementation.
  hasSpecificFamily_ = false;
  const String measureType(measure_.getImplementation()->getClassName());
  // Legendre factory
  if (measureType == "Uniform")
  {
    hasSpecificFamily_ = true;
    const Point parameter(measure_.getParameter());
    specificFamily_ = LegendreFactory();
    linear_ = 2.0 / (parameter[1] - parameter[0]);
    constant_ = -(parameter[0] + parameter[1]) / (parameter[1] - parameter[0]);
  }
  // Hermite factory
  if (measureType == "Normal")
  {
    hasSpecificFamily_ = true;
    const Point parameter(measure_.getParameter());
    specificFamily_ = HermiteFactory();
    linear_ = parameter[1];
    constant_ = parameter[0];
  }
  // HistogramPolynomial factory
  if (measureType == "Histogram")
  {
    hasSpecificFamily_ = true;
    const Point parameter(measure_.getParameter());
    const UnsignedInteger size = (parameter.getSize() - 1) / 2;
    const Scalar first = parameter[0];
    Point width(size);
    Point height(size);
    for (UnsignedInteger i = 0; i < size; ++i)
    {
      width[i] = parameter[2 * i + 1];
      height[i] = parameter[2 * i + 2];
    }
    specificFamily_ = HistogramPolynomialFactory(first, width, height);
  }
  // Chebychev factory
  if (measureType == "Arcsine")
  {
    hasSpecificFamily_ = true;
    specificFamily_ = ChebychevFactory();
  }
  // Jacobi (or Chebychev as a special case) factory
  if (measureType == "Beta")
  {
    hasSpecificFamily_ = true;
    const Point parameter(measure_.getParameter());
    const Scalar alpha = parameter[1] - parameter[0] - 1.0;
    const Scalar beta = parameter[0] - 1.0;
    const Scalar a = parameter[2];
    const Scalar b = parameter[3];
    linear_ = 2.0 / (b - a);
    constant_ = -(a + b) / (b - a);
    // Here we set directly the specific family as the reference distribution
    // of the family has a different type (Arcsine) than the given distribution
    if (alpha == -0.5 && beta == -0.5)
      specificFamily_ = ChebychevFactory();
    // Here we set directly the specific family as the reference distribution
    // of the family has a different type (Uniform) than the given distribution
    else if (alpha == 0.0 && beta == 0.0)
      specificFamily_ = LegendreFactory();
    else
      specificFamily_ = JacobiFactory(alpha, beta);
  }
  // Laguerre factory
  if (measureType == "Gamma")
  {
    hasSpecificFamily_ = true;
    const Point parameter(measure_.getParameter());
    specificFamily_ = LaguerreFactory(parameter[0] - 1.0);
    linear_ = 1.0 / parameter[1];
    constant_ = parameter[2];
  }
  if (measureType == "Exponential")
  {
    hasSpecificFamily_ = true;
    const Point parameter(measure_.getParameter());
    specificFamily_ = LaguerreFactory(0.0);
    linear_ = 1.0 / parameter[0];
    constant_ = parameter[1];
  }
  // Charlier factory
  if (measureType == "Poisson")
  {
    hasSpecificFamily_ = true;
    const Point parameter(measure_.getParameter());
    specificFamily_ = CharlierFactory(parameter[0]);
  }
  // Krawtchouk factory
  if (measureType == "Binomial")
  {
    hasSpecificFamily_ = true;
    const Point parameter(measure_.getParameter());
    specificFamily_ = KrawtchoukFactory(static_cast<UnsignedInteger>(parameter[0]), parameter[1]);
  }
  // Meixner factory
  if (measureType == "NegativeBinomial")
  {
    hasSpecificFamily_ = true;
    const Point parameter(measure_.getParameter());
    specificFamily_ = MeixnerFactory(parameter[0], parameter[1]);
  }
  if (!hasSpecificFamily_)
  {
    linear_ = measure_.getStandardDeviation()[0];
    constant_ = measure_.getMean()[0];
    measure_ = measure_ * (1.0 / linear_) + (-constant_ / linear_);
  }
}

/* String converter */
String StandardDistributionPolynomialFactory::__repr__() const
{
  OSS oss;
  oss << "class=" << getClassName()
      << " hasSpecificFamily=" << std::boolalpha << hasSpecificFamily_;
  if (hasSpecificFamily_)
  {
    oss << " specificFamily=" << specificFamily_
        << " linear=" << linear_
        << " constant=" << constant_;
  }
  else oss << " orthonormalization algorithm=" << orthonormalizationAlgorithm_;
  return oss;
}

/* Method save() stores the object through the StorageManager */
void StandardDistributionPolynomialFactory::save(Advocate & adv) const
{
  OrthogonalUniVariatePolynomialFactory::save(adv);
  adv.saveAttribute( "orthonormalizationAlgorithm_", orthonormalizationAlgorithm_ );
  adv.saveAttribute( "specificFamily_", specificFamily_ );
  adv.saveAttribute( "hasSpecificFamily_", hasSpecificFamily_ );
  adv.saveAttribute( "hasSpecificFamily_", hasSpecificFamily_ );
  adv.saveAttribute( "linear_", linear_ );
  adv.saveAttribute( "constant_", constant_ );
}


/* Method load() reloads the object from the StorageManager */
void StandardDistributionPolynomialFactory::load(Advocate & adv)
{
  OrthogonalUniVariatePolynomialFactory::load(adv);
  adv.loadAttribute( "orthonormalizationAlgorithm_", orthonormalizationAlgorithm_ );
  adv.loadAttribute( "specificFamily_", specificFamily_ );
  adv.loadAttribute( "hasSpecificFamily_", hasSpecificFamily_ );
  adv.loadAttribute( "linear_", linear_ );
  adv.loadAttribute( "constant_", constant_ );
}


END_NAMESPACE_OPENTURNS
