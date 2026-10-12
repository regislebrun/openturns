#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott

ot.TESTPREAMBLE()


# Instantiate one distribution object
distribution = ot.DiscreteCompoundDistribution(ot.Bernoulli(0.5), ot.Poisson(20.0))
upper_bound = int(distribution.getRange().getUpperBound()[0])
assert upper_bound == 62

# Compare with mathematically equal distribution
poisson_distribution = ot.Poisson(10.0)

for i in range(upper_bound):
    ott.assert_almost_equal(
        distribution.computePDF([i]), poisson_distribution.computePDF([i])
    )

assert distribution.getClassName() == "DiscreteCompoundDistribution"

# Is this distribution elliptical ?
assert not distribution.isElliptical()

# Is this distribution continuous ?
assert not distribution.isContinuous()

# Test for realization of distribution
oneRealization = distribution.getRealization()
ott.assert_almost_equal(oneRealization, [4.0], 0.0, 0.0)

# Test for sampling
size = 10000
oneSample = distribution.getSample(size)
ott.assert_almost_equal(oneSample[0], [12.0], 0.0, 0.0)
ott.assert_almost_equal(oneSample[size - 1], [18.0], 0.0, 0.0)
ott.assert_almost_equal(oneSample.computeMean(), [9.9932], 1e-5, 0.0)
ott.assert_almost_equal(
    oneSample.computeCovariance(), ot.CovarianceMatrix(1, [9.88734]), 1e-5, 0.0
)

# quantile
quantile = distribution.computeQuantile(0.95)
ott.assert_almost_equal(quantile, [15.0], 0.0, 0.0)
ott.assert_almost_equal(distribution.computeCDF(quantile), 0.951260, 1e-6, 1e-6)
ott.assert_almost_equal(distribution.computeEntropy(), 2.561410, 1e-6, 1e-6)
mean = distribution.getMean()
ott.assert_almost_equal(mean, [10.0], 1e-5, 0.0)
standardDeviation = distribution.getStandardDeviation()
ott.assert_almost_equal(standardDeviation, [3.16228], 1e-5, 0.0)
skewness = distribution.getSkewness()
ott.assert_almost_equal(skewness, [0.316228], 1e-5, 0.0)
kurtosis = distribution.getKurtosis()
ott.assert_almost_equal(kurtosis, [3.1], 1e-5, 0.0)
covariance = distribution.getCovariance()
ott.assert_almost_equal(covariance, ot.CovarianceMatrix(1, [10.0]), 1e-5, 0.0)
assert distribution.getStandardRepresentative().getImplementation().getClassName() == "DiscreteCompoundDistribution"

ot.Log.Show(ot.Log.TRACE)
validation = ott.DistributionValidation(distribution)
validation.skipParameters()  # see FiniteDiscreteDistribution
validation.skipInverseSurvival()  # numerical artifact ?
validation.run()

# Oversized range products are rejected before conversion
with ott.assert_raises(TypeError):
    ot.DiscreteCompoundDistribution(ot.Binomial(2**40), ot.Binomial(2**40))
