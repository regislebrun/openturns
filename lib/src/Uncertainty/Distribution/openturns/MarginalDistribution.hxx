//                                               -*- C++ -*-
/**
 *  @brief The class that implements a generic mechanism to extract marginal distributions
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
#ifndef OPENTURNS_MARGINALDISTRIBUTION_HXX
#define OPENTURNS_MARGINALDISTRIBUTION_HXX

#include "openturns/Distribution.hxx"
#include "openturns/DistributionImplementation.hxx"
#include "openturns/IntegrationAlgorithm.hxx"
#include "openturns/IteratedQuadrature.hxx"
#include "openturns/GaussKronrod.hxx"

BEGIN_NAMESPACE_OPENTURNS

/**
 * @class MarginalDistribution
 *
 * The class describes the probabilistic concept of marginal distributions
 */
class OT_API MarginalDistribution
  : public DistributionImplementation
{
  CLASSNAME
public:

  /** Default constructor */
  MarginalDistribution();

  /** Parameters constructor */
  MarginalDistribution(const Distribution & distribution,
                       const Indices & indices);

  /** Comparison operator */
  using DistributionImplementation::operator ==;
  Bool operator ==(const MarginalDistribution & other) const;
protected:
  Bool equals(const DistributionImplementation & other) const override;
public:

  /** String converter */
  String __repr__() const override;
  String __str__(const String & offset = "") const override;


  /** Distribution accessor */
  void setDistribution(const Distribution & distribution);
  Distribution getDistribution() const;

  /** Indices accessor */
  void setIndices(const Indices & indices);
  Indices getIndices() const;

  /** UsePDF accessor */
  void setUsePDF(const Bool usePDF);
  Bool getUsePDF() const;

  /** Integration algorithm accessor */
  void setIntegrationAlgorithm(const IntegrationAlgorithm & algo);
  IntegrationAlgorithm getIntegrationAlgorithm() const;

  /** Get the support of a discrete distribution that intersect a given interval */
  using DistributionImplementation::getSupport;
  Sample getSupport() const override;
  Sample getSupport(const Interval & interval) const override;

  /** Get the discrete probability levels */
  using DistributionImplementation::getProbabilities;
  Point getProbabilities() const override;

private:
  /** Set the distribution and the indices in one shot */
  void setDistributionAndIndices(const Distribution & distribution,
                                 const Indices & indices);
public:

  /* Here is the interface that all derived class must implement */

  /** Virtual constructor */
  MarginalDistribution * clone() const override;

  /** Get one realization of the MarginalDistribution */
  Point getRealization() const override;
  Sample getSample(const UnsignedInteger size) const override;

  /** Get the PDF of the MarginalDistribution */
  using DistributionImplementation::computePDF;
  Scalar computePDF(const Point & point) const override;

  /** Get the CDF of the MarginalDistribution */
  using DistributionImplementation::computeCDF;
  Scalar computeCDF(const Point & point) const override;

  /** Get the survival function of the MarginalDistribution */
  using DistributionImplementation::computeSurvivalFunction;
  Scalar computeSurvivalFunction(const Point & point) const override;

  /** Get the probability content of an interval */
  Scalar computeProbability(const Interval & interval) const override;

  /** Get the standard deviation of the distribution */
  Point getStandardDeviation() const override;

  /** Get the skewness of the distribution */
  Point getSkewness() const override;

  /** Get the kurtosis of the distribution */
  Point getKurtosis() const override;

  /** Get the Spearman correlation of the distribution */
  CorrelationMatrix getSpearmanCorrelation() const override;

  /** Get the Kendall concordance of the distribution */
  CorrelationMatrix getKendallTau() const override;

  /** Get the i-th marginal distribution */
  Distribution getMarginal(const UnsignedInteger i) const override;

  /** Get the distribution of the marginal distribution corresponding to indices dimensions */
  Distribution getMarginal(const Indices & indices) const override;

  /** Get the standard distribution */
  Distribution getStandardDistribution() const override;

  /** Tell if the distribution has independent copula */
  Bool hasIndependentCopula() const override;

  /** Parameters value and description accessor */
  PointWithDescriptionCollection getParametersCollection() const override;

  /** Tell if the distribution has elliptical copula */
  Bool hasEllipticalCopula() const override;

  /** Check if the distribution is elliptical */
  Bool isElliptical() const override;

  /** Check if the distribution is continuous */
  Bool isContinuous() const override;

  /** Check if the distribution is discrete */
  Bool isDiscrete() const override;

  /** Tell if the distribution is integer valued */
  Bool isIntegral() const override;

  /** Parameters accessors */
  void setParameter(const Point & parameter) override;
  Point getParameter() const override;
  Description getParameterDescription() const override;

  /** Method save() stores the object through the StorageManager */
  void save(Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(Advocate & adv) override;


protected:
  /** Compute the mean of the distribution */
  void computeMean() const override;

  /** Compute the covariance of the distribution */
  void computeCovariance() const override;

  /** Expand the given marginal point to the underlying distribution argument point */
  Point expandPoint(const Point & point,
                    const Bool upper = true) const;

  /** Reduce the given point to the marginal point */
  Point reducePoint(const Point & point) const;

  /** The distribution under the MarginalDistribution */
  Distribution distribution_;

  /** The indices of the MarginalDistribution */
  Indices indices_;

  /** The lower bound of the underlying distribution */
  Point lowerBound_;

  /** The upper bound of the underlying distribution */
  Point upperBound_;

  /** Flag to tell if the computations have to be done using the CDF or the PDF
      of the underlying distribution */
  Bool usePDF_ = ResourceMap::GetAsBool("MarginalDistribution-UsePDF");

  /** Integration algorithm used to compute the PDF. By default it is an IteratedQuadrature. */
  IntegrationAlgorithm integrationAlgorithm_ = IteratedQuadrature(
        GaussKronrod(ResourceMap::GetAsUnsignedInteger("MarginalDistribution-MaximumSubIntervals"),
                     ResourceMap::GetAsScalar("MarginalDistribution-MaximumError"),
                     GaussKronrod::GetRuleFromName(ResourceMap::GetAsString("MarginalDistribution-Rule"))
                    ));

}; /* class MarginalDistribution */


END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_MARGINALDISTRIBUTION_HXX */
