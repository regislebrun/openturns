#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott

ot.TESTPREAMBLE()

ot.ResourceMap.SetAsUnsignedInteger("MeixnerDistribution-CDFDiscretization", 100)
ot.ResourceMap.SetAsScalar("MeixnerDistribution-MaximumAbsoluteError", 1.0e-6)
ot.ResourceMap.SetAsScalar("MeixnerDistribution-MaximumRelativeError", 1.0e-6)
ot.ResourceMap.SetAsScalar("MeixnerDistribution-MaximumConstraintError", 1.0e-6)
ot.ResourceMap.SetAsScalar("MeixnerDistribution-MaximumObjectiveError", 1.0e-6)
distribution = ot.MeixnerDistribution(1.5, 0.5, 2.5, -0.5)
assert distribution.getClassName() == "MeixnerDistribution"
size = 1000
sample = distribution.getSample(size)
factory = ot.MeixnerDistributionFactory()
covariance = ot.CovarianceMatrix()
# Distribution estimatedDistribution(factory.build(sample, covariance))
estimatedDistribution = factory.build(sample)
ott.assert_almost_equal(estimatedDistribution.getBeta(), 1.77623, 1e-5, 0.0)
ott.assert_almost_equal(estimatedDistribution.getAlpha(), 0.316056, 1e-5, 0.0)
ott.assert_almost_equal(estimatedDistribution.getDelta(), 1.77799, 1e-5, 0.0)
ott.assert_almost_equal(estimatedDistribution.getGamma(), -0.0426543, 1e-5, 0.0)
estimatedDistribution = factory.build()
ott.assert_almost_equal(estimatedDistribution.getBeta(), 1.0, 0.0, 0.0)
ott.assert_almost_equal(estimatedDistribution.getAlpha(), 0.0, 0.0, 0.0)
ott.assert_almost_equal(estimatedDistribution.getDelta(), 1.0, 0.0, 0.0)
ott.assert_almost_equal(estimatedDistribution.getGamma(), 0.0, 0.0, 0.0)
estimatedDistribution = factory.build(distribution.getParameter())
ott.assert_almost_equal(estimatedDistribution.getParameter(), distribution.getParameter(), 0.0, 0.0)
estimatedMeixnerDistribution = factory.buildAsMeixnerDistribution(sample)
ott.assert_almost_equal(estimatedMeixnerDistribution.getBeta(), 1.77623, 1e-5, 0.0)
ott.assert_almost_equal(estimatedMeixnerDistribution.getAlpha(), 0.316056, 1e-5, 0.0)
ott.assert_almost_equal(estimatedMeixnerDistribution.getDelta(), 1.77799, 1e-5, 0.0)
ott.assert_almost_equal(estimatedMeixnerDistribution.getGamma(), -0.0426543, 1e-5, 0.0)
estimatedMeixnerDistribution = factory.buildAsMeixnerDistribution()
ott.assert_almost_equal(estimatedMeixnerDistribution.getBeta(), 1.0, 0.0, 0.0)
estimatedMeixnerDistribution = factory.buildAsMeixnerDistribution(
    distribution.getParameter()
)
ott.assert_almost_equal(estimatedMeixnerDistribution.getParameter(), distribution.getParameter(), 0.0, 0.0)

# Moment matching preserves the sample variance and mean by construction
momentEstimate = factory.buildAsMeixnerDistribution(sample)
ott.assert_almost_equal(
    momentEstimate.getCovariance()[0, 0], sample.computeVariance()[0], 1e-9, 0.0
)
ott.assert_almost_equal(
    momentEstimate.getMean()[0], sample.computeMean()[0], 1e-9, 0.0
)
