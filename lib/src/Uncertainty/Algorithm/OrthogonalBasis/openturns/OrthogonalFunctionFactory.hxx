//                                               -*- C++ -*-
/**
 *  @brief This is an abstract class for orthogonal basis
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
#ifndef OPENTURNS_ORTHOGONALFUNCTIONFACTORY_HXX
#define OPENTURNS_ORTHOGONALFUNCTIONFACTORY_HXX

#include "openturns/BasisImplementation.hxx"
#include "openturns/Function.hxx"
#include "openturns/Distribution.hxx"
#include "openturns/EnumerateFunction.hxx"

BEGIN_NAMESPACE_OPENTURNS

class OrthogonalBasis;

/**
 * @class OrthogonalFunctionFactory
 *
 * This is an abstract class for orthogonal basis
 */

class OT_API OrthogonalFunctionFactory
  : public BasisImplementation
{
  CLASSNAME
public:


  /** Default constructor */
  OrthogonalFunctionFactory();

  /** Parameter constructor */
  OrthogonalFunctionFactory(const Distribution & measure);

  /** Build the Function of the given index */
  Function build(const UnsignedInteger index) const override;

  /** Build the Function of the given multi-indices */
  virtual Function build(const Indices & indices) const;

  /** Return the measure upon which the basis is orthogonal */
  virtual Distribution getMeasure() const;

  /** Return the enumerate function that translate unidimensional indices into multidimensional indices */
  virtual EnumerateFunction getEnumerateFunction() const;

  /** Virtual constructor */
  OrthogonalFunctionFactory * clone() const override;

  /** Get the function factory corresponding to marginal input indices */
  virtual OrthogonalBasis getMarginal(const Indices & indices) const;

  Bool isOrthogonal() const override;

  Bool isTensorProduct() const override;

  /** String converter */
  String __repr__() const override;

  /** Method save() stores the object through the StorageManager */
  void save(Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(Advocate & adv) override;


protected:
  /** The measure that defines the scalar product */
  Distribution measure_;

private:

} ; /* class OrthogonalFunctionFactory */


END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_ORTHOGONALFUNCTIONFACTORY_HXX */
