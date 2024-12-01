//                                               -*- C++ -*-
/**
 *  @brief Abstract top-level view of an ratioOfUniformsExperiment plane
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
#include "openturns/RatioOfUniformsExperiment.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/RandomGenerator.hxx"
#include "openturns/SobolSequence.hxx"
#include "openturns/OptimizationAlgorithm.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(RatioOfUniformsExperiment)

static const Factory<RatioOfUniformsExperiment> Factory_RatioOfUniformsExperiment;

/* Default constructor */
RatioOfUniformsExperiment::RatioOfUniformsExperiment()
  : WeightedExperimentImplementation()
{
  // Nothing to do
}

/* Constructor with parameters */
RatioOfUniformsExperiment::RatioOfUniformsExperiment(const UnsignedInteger size)
  : WeightedExperimentImplementation(size)
{
  // Nothing to do
}

/* Constructor with parameters */
RatioOfUniformsExperiment::RatioOfUniformsExperiment(const Distribution & distribution,
    const UnsignedInteger size)
  : WeightedExperimentImplementation(distribution, size)
{
  // Prepare the ratio of uniforms method
  candidateNumber_ = ResourceMap::GetAsUnsignedInteger("RatioOfUniformsExperiment-CandidateNumber");
  optimizationAlgorithm_ = OptimizationAlgorithm::GetByName(ResourceMap::GetAsString("RatioOfUniformsExperiment-OptimizationAlgorithm"));
  setDistribution(distribution_);
}

/* Virtual constructor */
RatioOfUniformsExperiment * RatioOfUniformsExperiment::clone() const
{
  return new RatioOfUniformsExperiment(*this);
}

/* String converter */
String RatioOfUniformsExperiment::__repr__() const
{
  OSS oss;
  oss << "class=" << GetClassName()
      << " name=" << getName()
      << " distribution=" << distribution_
      << " size=" << size_;
  return oss;
}

/* Optimization algorithm accessor */
void RatioOfUniformsExperiment::setOptimizationAlgorithm(const OptimizationAlgorithm & optimizationAlgorithm)
{
  optimizationAlgorithm_ = optimizationAlgorithm;
}

OptimizationAlgorithm RatioOfUniformsExperiment::getOptimizationAlgorithm() const
{
  return optimizationAlgorithm_;
}

/* Candidate number accessor */
void RatioOfUniformsExperiment::setCandidateNumber(const UnsignedInteger candidateNumber)
{
  if (candidateNumber == 0)
    throw InvalidArgumentException(HERE) << "Error: the candidate number must be at least 1";
  candidateNumber_ = candidateNumber;
}

UnsignedInteger RatioOfUniformsExperiment::getCandidateNumber() const
{
  return candidateNumber_;
}

class RatioOfUniformsExperimentUBoundEvaluation : public EvaluationImplementation
{
public:
  RatioOfUniformsExperimentUBoundEvaluation(const Distribution & distribution, const Scalar r)
    : EvaluationImplementation()
    , distribution_(distribution)
    , r_(r)
  {
    // Nothing to do
  }

  RatioOfUniformsExperimentUBoundEvaluation * clone() const override
  {
    return new RatioOfUniformsExperimentUBoundEvaluation(*this);
  }

  UnsignedInteger getInputDimension() const override
  {
    return distribution_.getDimension();
  }

  UnsignedInteger getOutputDimension() const override
  {
    return 1;
  }

  Point operator()(const Point & inP) const override
  {
    const UnsignedInteger dimension = distribution_.getDimension();
    Scalar result = distribution_.computeLogPDF(inP) / (1.0 + r_ * dimension);
    result = std::max(-SpecFunc::LogMaxScalar, result);
    return {result};
  }

private:
  Distribution distribution_;
  Scalar r_ = 0.0;

};

class RatioOfUniformsExperimentVBoundEvaluation : public EvaluationImplementation
{
public:
  RatioOfUniformsExperimentVBoundEvaluation(const Distribution & distribution, const Scalar r)
    : EvaluationImplementation()
    , distribution_(distribution)
    , r_(r)
  {
    // Nothing to do
  }

  RatioOfUniformsExperimentVBoundEvaluation * clone() const override
  {
    return new RatioOfUniformsExperimentVBoundEvaluation(*this);
  }

  UnsignedInteger getInputDimension() const override
  {
    return distribution_.getDimension();
  }

  UnsignedInteger getOutputDimension() const override
  {
    return distribution_.getDimension();
  }

  Point operator()(const Point & inP) const override
  {
    const UnsignedInteger dimension = distribution_.getDimension();
    const Scalar value = distribution_.computeLogPDF(inP) * r_ / (1.0 + r_ * dimension);
    Point result(dimension, value);
    for (UnsignedInteger i = 0; i < dimension; ++ i)
    {
      result[i] += std::log(std::abs(inP[i]));
      result[i] = std::max(-SpecFunc::LogMaxScalar, result[i]);
    }
    return result;
  }

private:
  Distribution distribution_;
  Scalar r_ = 0.0;
};

/* Distribution accessor */
void RatioOfUniformsExperiment::setDistribution(const Distribution & distribution)
{
  if (!distribution.isContinuous())
    throw InvalidArgumentException(HERE) << "Error: the ratio of uniforms algorithm works only with continuous distributions, here distribution=" << distribution;
  // r_ is a free parameter, could be optimized to maximize the acceptance ratio
  const UnsignedInteger dimension = distribution.getDimension();
  const Interval bounds(distribution.getRange());
  const Point lb(bounds.getLowerBound());
  const Point ub(bounds.getUpperBound());

  // find a feasible starting point
  SobolSequence sequence(dimension);
  Point start;
  for (UnsignedInteger k = 0; k < candidateNumber_; ++ k)
    {
      Point candidate(sequence.generate());
      for (UnsignedInteger j = 0; j < dimension; ++ j)
        candidate[j] = lb[j] + candidate[j] * (ub[j] - lb[j]);
      if (SpecFunc::IsNormal(distribution.computeLogPDF(candidate)))
	{
	  start = candidate;
	  break;
	}
    } // for k
  if (!start.getDimension())
    throw InternalException(HERE) << "Could not find a feasible starting point to initialize ration of uniforms U sup";
  
  // First, the upper bound on U
  const Function objectiveU(new RatioOfUniformsExperimentUBoundEvaluation(distribution, r_));
  OptimizationProblem problemU(objectiveU);
  problemU.setMinimization(false);
  problemU.setBounds(bounds);
  optimizationAlgorithm_.setProblem(problemU);
  optimizationAlgorithm_.setStartingPoint(start);
  optimizationAlgorithm_.run();
  supU_ = std::exp(optimizationAlgorithm_.getResult().getOptimalValue()[0]);
  LOGDEBUG(OSS() << "supU_=" << supU_ << " u*=" << optimizationAlgorithm_.getResult().getOptimalPoint());
  
  // Second, the lower and upper bounds on V
  const Function objectiveV(new RatioOfUniformsExperimentVBoundEvaluation(distribution, r_));
  infV_.resize(dimension);
  supV_.resize(dimension);
  const Point zero(dimension, 0.0);
  for (UnsignedInteger i = 0; i < dimension; ++ i)
    {
      const Function objectiveVI(objectiveV.getMarginal(i));
      OptimizationProblem problemVI(objectiveVI);
      problemVI.setMinimization(false);
      if (ub[i] > 0.0)
	{
	  // find a feasible starting point in [0, ub]
	  start.clear();
	  for (UnsignedInteger k = 0; k < candidateNumber_; ++ k)
	    {
	      Point candidate(sequence.generate());
	      for (UnsignedInteger j = 0; j < dimension; ++ j)
		candidate[j] = candidate[j] * ub[j];
	      if (SpecFunc::IsNormal(distribution.computeLogPDF(candidate)))
		{
		  start = candidate;
		  break;
		}
	    } // for k
	  if (!start.getDimension())
	    throw InternalException(HERE) << "Could not find a feasible starting point to initialize ration of uniforms V sup";
	  problemVI.setBounds(Interval(zero, ub));
	  optimizationAlgorithm_.setProblem(problemVI);
	  optimizationAlgorithm_.setStartingPoint(start);
	  optimizationAlgorithm_.run();
	  supV_[i] = std::exp(optimizationAlgorithm_.getResult().getOptimalValue()[0]);
	  LOGDEBUG(OSS() << "supV_[" << i << "]=" << supV_[i] << " v*=" << optimizationAlgorithm_.getResult().getOptimalPoint());
	} // if ub[i] > 0.0
      if (lb[i] < 0.0)
      {
        // find a feasible starting point in [lb, 0]
        start.clear();
        for (UnsignedInteger k = 0; k < candidateNumber_; ++ k)
	  {
	    Point candidate(sequence.generate());
	    for (UnsignedInteger j = 0; j < dimension; ++ j)
	      candidate[j] = candidate[j] * lb[j];
	    if (SpecFunc::IsNormal(distribution.computeLogPDF(candidate)))
	      {
		start = candidate;
		break;
	      }
	  } // for k
        if (!start.getDimension())
          throw InternalException(HERE) << "Could not find a feasible starting point to initialize ration of uniforms V inf";
        problemVI.setBounds(Interval(lb, zero));
        optimizationAlgorithm_.setProblem(problemVI);
        optimizationAlgorithm_.setStartingPoint(start);
        optimizationAlgorithm_.run();
        infV_[i] = -std::exp(optimizationAlgorithm_.getResult().getOptimalValue()[0]);
        LOGDEBUG(OSS() << "infV_[" << i << "]=" << infV_[i] << " v*=" << optimizationAlgorithm_.getResult().getOptimalPoint());
      } // if lb[i] < 0.0
    } // for i
  WeightedExperimentImplementation::setDistribution(distribution);
}

/* Sample generation */
Sample RatioOfUniformsExperiment::generateWithWeights(Point & weights) const
{
  weights = Point(size_, 1.0 / size_);
  if (!infV_.getSize())
    throw InvalidArgumentException(HERE) << "RatioOfUniformsExperiment was not initialized. Call setDistribution() to fix it.";

  // Now, the sampling using rejection
  const UnsignedInteger dimension = distribution_.getDimension();
  Sample sample(size_, dimension);
  Point result(dimension);
  for (UnsignedInteger n = 0; n < size_; ++n)
    {
      Bool accepted = false;
      while (!accepted)
	{
	  const Scalar u = supU_ * RandomGenerator::Generate();
	    const Scalar ur = std::pow(u, r_);
	    for (UnsignedInteger i = 0; i < dimension; ++ i)
	      result[i] = (infV_[i] + (supV_[i] - infV_[i]) * RandomGenerator::Generate()) / ur;
	    accepted = (1.0 + r_ * dimension) * std::log(u) <= distribution_.computeLogPDF(result);
	  } // !accepted
      sample[n] = result;
    } // for n
  return sample;
}

END_NAMESPACE_OPENTURNS
