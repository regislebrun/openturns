//                                               -*- C++ -*-
/**
 *  @brief Abstract top-level view of a RatioOfUniformsCarloExperiment generator
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
#ifndef OPENTURNS_RATIOOFUNIFORMSEXPERIMENT_HXX
#define OPENTURNS_RATIOOFUNIFORMSEXPERIMENT_HXX

#include "openturns/WeightedExperimentImplementation.hxx"
#include "openturns/OptimizationAlgorithm.hxx"

BEGIN_NAMESPACE_OPENTURNS



/**
 * @class RatioOfUniformsExperiment
 *
 * The class describes the probabilistic concept of ratioOfUniformsExperiment generator
 */
class OT_API RatioOfUniformsExperiment
  : public WeightedExperimentImplementation
{
  CLASSNAME
public:


  /** Default constructor */
  RatioOfUniformsExperiment();

  /** Parameters constructor */
  explicit RatioOfUniformsExperiment(const UnsignedInteger size);

  /** Parameters constructor */
  RatioOfUniformsExperiment(const Distribution & distribution,
                       const UnsignedInteger size);

  /** Virtual constructor */
  RatioOfUniformsExperiment * clone() const override;

  /** String converter */
  String __repr__() const override;

  /** Distribution accessor */
  void setDistribution(const Distribution & distribution) override;

  /** Optimization algorithm accessor */
  void setOptimizationAlgorithm(const OptimizationAlgorithm & optimizationAlgorithm);
  OptimizationAlgorithm getOptimizationAlgorithm() const;

  /** Candidate number accessor */
  void setCandidateNumber(const UnsignedInteger candidateNumber);
  UnsignedInteger getCandidateNumber() const;

  /* Here is the interface that all derived class must implement */

  /** Sample generation */
  Sample generateWithWeights(Point & weightsOut) const override;

protected:

private:
  Scalar r_ = 1.0;
  Scalar supU_ = 0.0;
  Point infV_;
  Point supV_;
  OptimizationAlgorithm optimizationAlgorithm_;
  UnsignedInteger candidateNumber_ = 1;
}; /* class RatioOfUniformsExperiment */


END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_RATIOOFUNIFORMSEXPERIMENT_HXX */
