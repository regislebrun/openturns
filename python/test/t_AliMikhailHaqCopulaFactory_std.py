#! /usr/bin/env python

import math

import openturns as ot
import openturns.testing as ott

ot.TESTPREAMBLE()


distribution = ot.AliMikhailHaqCopula(0.5)
size = 1000
sample = distribution.getSample(size)
factory = ot.AliMikhailHaqCopulaFactory()
estimatedDistribution = factory.build(sample)
assert distribution.getClassName() == "AliMikhailHaqCopula"
ott.assert_almost_equal(distribution.getTheta(), 0.5, 0.0, 0.0)
ott.assert_almost_equal(estimatedDistribution.getTheta(), 0.433481, 1e-5, 0.0)
estimatedDistribution = factory.build()
ott.assert_almost_equal(estimatedDistribution.getTheta(), 0.5, 0.0, 0.0)
estimatedDistribution = factory.build(distribution.getParameter())
ott.assert_almost_equal(estimatedDistribution.getTheta(), 0.5, 0.0, 0.0)
estimatedAliMikhailHaqCopula = factory.buildAsAliMikhailHaqCopula(sample)
ott.assert_almost_equal(estimatedAliMikhailHaqCopula.getTheta(), 0.433481, 1e-5, 0.0)
estimatedAliMikhailHaqCopula = factory.buildAsAliMikhailHaqCopula()
ott.assert_almost_equal(estimatedAliMikhailHaqCopula.getTheta(), 0.5, 0.0, 0.0)
estimatedAliMikhailHaqCopula = factory.buildAsAliMikhailHaqCopula(
    distribution.getParameter()
)
ott.assert_almost_equal(estimatedAliMikhailHaqCopula.getTheta(), 0.5, 0.0, 0.0)

# A non-finite bisection tolerance falls back to the default precision
ot.RandomGenerator.SetSeed(0)
toleranceSample = ot.AliMikhailHaqCopula(0.5).getSample(200)
defaultEstimated = factory.buildAsAliMikhailHaqCopula(toleranceSample)
ot.ResourceMap.SetAsScalar("AliMikhailHaqCopulaFactory-ThetaEpsilon", float("inf"))
try:
    fallbackEstimated = factory.buildAsAliMikhailHaqCopula(toleranceSample)
    theta = fallbackEstimated.getTheta()
    assert math.isfinite(theta)
    assert -1.0 <= theta <= 1.0
    # without fallback the unbracketed edge would be returned instead
    ott.assert_almost_equal(
        theta, defaultEstimated.getTheta(), 1e-9, 0.0
    )
finally:
    ot.ResourceMap.SetAsScalar("AliMikhailHaqCopulaFactory-ThetaEpsilon", 1.0e-14)
