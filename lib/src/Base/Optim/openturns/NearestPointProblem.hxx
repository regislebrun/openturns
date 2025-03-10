//                                               -*- C++ -*-
/**
 *  @brief NearestPointProblem allows one to describe an optimization problem
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
#ifndef OPENTURNS_NEARESTPOINTPROBLEM_HXX
#define OPENTURNS_NEARESTPOINTPROBLEM_HXX

#include "openturns/OptimizationProblemImplementation.hxx"

BEGIN_NAMESPACE_OPENTURNS

/**
 * @class NearestPointProblem
 * NearestPointProblem allows one to describe an optimization problem
 */

class OT_API NearestPointProblem
  : public OptimizationProblemImplementation
{

  CLASSNAME

public:

  /** Default constructor */
  NearestPointProblem();

  /** Constructor with parameters */
  NearestPointProblem(const Function & levelFunction,
                      Scalar levelValue);

  /** Virtual constructor */
  NearestPointProblem * clone() const override;

  /** Level function accessor */
  Function getLevelFunction() const override;
  void setLevelFunction(const Function & levelFunction) override;
  Bool hasLevelFunction() const override;

  /** Level value accessor */
  Scalar getLevelValue() const override;
  void setLevelValue(Scalar levelValue) override;

  /** String converter */
  String __repr__() const override;

  /** Method save() stores the object through the StorageManager */
  void save(Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(Advocate & adv) override;

private:
  void clearLevelFunction();
  void setNearestPointConstraints();

  // The level function, for nearest point problems
  Function levelFunction_;

  // The level value, for nearest point problems
  Scalar levelValue_;

} ; /* class NearestPointProblem */


END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_NEARESTPOINTPROBLEM_HXX */
