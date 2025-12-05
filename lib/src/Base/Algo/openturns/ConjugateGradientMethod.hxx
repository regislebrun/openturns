//                                               -*- C++ -*-
/**
 *  @brief ConjugateGradient iterative LS solver
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
#ifndef OPENTURNS_CONJUGATEGRADIENTMETHOD_HXX
#define OPENTURNS_CONJUGATEGRADIENTMETHOD_HXX

#include "openturns/LeastSquaresMethodImplementation.hxx"

BEGIN_NAMESPACE_OPENTURNS

/**
 * @class ConjugateGradientMethod
 *
 * ConjugateGradient decomposition based LS solver
 */
class OT_API ConjugateGradientMethod
  : public LeastSquaresMethodImplementation
{
  CLASSNAME
public:

  /** Default constructor */
  ConjugateGradientMethod();

  /** Parameters constructor */
  ConjugateGradientMethod(const DesignProxy & proxy,
                 const Point & weight,
                 const Indices & indices);

  /** Parameters constructor */
  ConjugateGradientMethod(const DesignProxy & proxy,
                 const Indices & indices);

  /** Parameters constructor */
  explicit ConjugateGradientMethod(const Matrix & matrix);

  /** Virtual constructor */
  ConjugateGradientMethod * clone() const override;

  /** String converter */
  String __repr__() const override;

  /** Starting point accessor */
  void setStartingPoint(const Point & startingPoint);
  Point getStartingPoint() const;

  /** Precision accessor */
  void setPrecision(const Scalar precision);
  Scalar getPrecision() const;

  /** Reinitialize the algorithm */
  void trashDecomposition() override;

  /** Solve least-squares problem, ie x=\argmin |Mx-b|^2 */
  Point solve(const Point & rhs) override;
  Point solveNormal(const Point & rhs) override;

  /** Method save() stores the object through the StorageManager */
  void save(Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(Advocate & adv) override;

protected:
  // Starting point for the conjugate gradient
  Point startingPoint_;

  // Maximum norm of the conjugate residual
  Scalar precision_;
}; /* class ConjugateGradientMethod */


END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_CONJUGATEGRADIENTMETHOD_HXX */
