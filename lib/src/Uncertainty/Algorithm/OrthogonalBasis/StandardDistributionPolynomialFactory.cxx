//                                               -*- C++
/**
 *  @brief StandardDistributionPolynomialFactory (deprecated wrapper)
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

#include "openturns/StandardDistributionPolynomialFactory.hxx"
#include "openturns/PersistentObjectFactory.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(StandardDistributionPolynomialFactory)

static const Factory<StandardDistributionPolynomialFactory> Factory_StandardDistributionPolynomialFactory;


StandardDistributionPolynomialFactory::StandardDistributionPolynomialFactory()
  : UniVariateDistributionPolynomialFactory()
{
  LOGWARN("StandardDistributionPolynomialFactory is deprecated");
}


StandardDistributionPolynomialFactory::StandardDistributionPolynomialFactory(const Distribution & distribution)
  : UniVariateDistributionPolynomialFactory(distribution)
{
  LOGWARN("StandardDistributionPolynomialFactory is deprecated");
}


StandardDistributionPolynomialFactory::StandardDistributionPolynomialFactory(const OrthonormalizationAlgorithm & orthonormalizationAlgorithm)
  : UniVariateDistributionPolynomialFactory(orthonormalizationAlgorithm)
{
  LOGWARN("StandardDistributionPolynomialFactory is deprecated");
}


StandardDistributionPolynomialFactory * StandardDistributionPolynomialFactory::clone() const
{
  return new StandardDistributionPolynomialFactory(*this);
}

/* Comparison operators */
Bool StandardDistributionPolynomialFactory::operator ==(const StandardDistributionPolynomialFactory & other) const
{
  if (this == &other) return true;
  if (!hasEqualBase(other)) return false;
  if (hasSpecificFamily_ != other.hasSpecificFamily_) return false;
  if (hasSpecificFamily_ && !(specificFamily_ == other.specificFamily_)) return false;
  return true;
}

Bool StandardDistributionPolynomialFactory::equals(const OrthogonalUniVariatePolynomialFactory & other) const
{
  const StandardDistributionPolynomialFactory * p_other = dynamic_cast<const StandardDistributionPolynomialFactory *>(&other);
  return p_other && (*this == *p_other);
}

END_NAMESPACE_OPENTURNS
