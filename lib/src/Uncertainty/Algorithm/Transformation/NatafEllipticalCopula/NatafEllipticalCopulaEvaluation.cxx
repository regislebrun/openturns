
//                                               -*- C++ -*-
/**
 *  @brief Class for the Nataf transformation evaluation for elliptical
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
#include "openturns/OTprivate.hxx"
#include "openturns/NatafEllipticalCopulaEvaluation.hxx"
#include "openturns/PersistentObjectFactory.hxx"

BEGIN_NAMESPACE_OPENTURNS

/*
 * @class NatafEllipticalCopulaEvaluation
 *
 * This class offers an interface for the Nataf function for elliptical copula
 */

CLASSNAMEINIT(NatafEllipticalCopulaEvaluation)

static const Factory<NatafEllipticalCopulaEvaluation> Factory_NatafEllipticalCopulaEvaluation;

/* Default constructor */
NatafEllipticalCopulaEvaluation::NatafEllipticalCopulaEvaluation()
  : EvaluationImplementation()
{
  // Nothing to do
}

/* Parameter constructor */
NatafEllipticalCopulaEvaluation::NatafEllipticalCopulaEvaluation(const Distribution & standardDistribution,
    const TriangularMatrix & inverseCholesky)
  : EvaluationImplementation()
  , standardDistribution_(standardDistribution)
  , inverseCholesky_(inverseCholesky)
{
  Description description(Description::BuildDefault(inverseCholesky.getDimension(), "x"));
  description.add(Description::BuildDefault(inverseCholesky.getDimension(), "y"));
  setDescription(description);
}

/* Virtual constructor */
NatafEllipticalCopulaEvaluation * NatafEllipticalCopulaEvaluation::clone() const
{
  return new NatafEllipticalCopulaEvaluation(*this);
}

/* String converter */
String NatafEllipticalCopulaEvaluation::__repr__() const
{
  OSS oss;
  oss << "class=" << NatafEllipticalCopulaEvaluation::GetClassName()
      << " description=" << getDescription()
      << " standardDistribution=" << standardDistribution_
      << " inverseCholesky=" << inverseCholesky_;

  return oss;
}

String NatafEllipticalCopulaEvaluation::__str__(const String & ) const
{
  OSS oss(false);
  oss << NatafEllipticalCopulaEvaluation::GetClassName()
      << "(Copula(inverseCholesky=" << inverseCholesky_ << ", E=" << standardDistribution_.getMarginal(0) << ")->" << standardDistribution_ << ")";

  return oss;
}

/*
 * Evaluation
 * This function transforms an elliptical copula with correlation matrix R into the associated standard elliptical distribution:
 * Yi(x) = Q(xi), where Q = F^{-1} and F is the CDF of the standard elliptical distribution
 * T(x) = G.Y(x), where G = L^{-1} and L is the Cholesky factor of R: L.L^t = R, L is lower triangular
 */
Point NatafEllipticalCopulaEvaluation::operator () (const Point & inP) const
{
  const UnsignedInteger dimension = getOutputDimension();
  Point result(dimension);
  const Distribution standardMarginal(standardDistribution_.getMarginal(0));
  // First, filter the common marginal distribution
  for (UnsignedInteger i = 0; i < dimension; ++ i)
    result[i] = standardMarginal.computeScalarQuantile(inP[i]);
  // Second, decorrelate the components
  result = inverseCholesky_ * result;
  callsNumber_.increment();
  return result;
}

/* Gradient according to the marginal parameters. Currently, the dependence parameters are not taken into account. */

Matrix NatafEllipticalCopulaEvaluation::parameterGradient(const Point & ) const
{
  return Matrix(0, getInputDimension());
}

/* Accessor for input point dimension */
UnsignedInteger NatafEllipticalCopulaEvaluation::getInputDimension() const
{
  return inverseCholesky_.getDimension();
}

/* Accessor for output point dimension */
UnsignedInteger NatafEllipticalCopulaEvaluation::getOutputDimension() const
{
  return inverseCholesky_.getDimension();
}

/* Method save() stores the object through the StorageManager */
void NatafEllipticalCopulaEvaluation::save(Advocate & adv) const
{
  EvaluationImplementation::save(adv);
  adv.saveAttribute( "standardDistribution_", standardDistribution_ );
  adv.saveAttribute( "inverseCholesky_", inverseCholesky_ );
}

/* Method load() reloads the object from the StorageManager */
void NatafEllipticalCopulaEvaluation::load(Advocate & adv)
{
  EvaluationImplementation::load(adv);
  adv.loadAttribute( "standardDistribution_", standardDistribution_ );
  if (adv.hasAttribute("inverseCholesky_"))
    adv.loadAttribute("inverseCholesky_", inverseCholesky_);
  else
    adv.loadAttribute("cholesky_", inverseCholesky_); // OT < 1.24
}

END_NAMESPACE_OPENTURNS
