#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import os

ot.TESTPREAMBLE()

# Instantiate one distribution object
distribution = ot.Binomial(15, 0.7)
assert distribution.getClassName() == "Binomial"

# Is this distribution elliptical ?
assert not distribution.isElliptical()

# Is this distribution continuous ?
assert not distribution.isContinuous()

# Test for realization of distribution
oneRealization = distribution.getRealization()
ott.assert_almost_equal(oneRealization, [10.0], 0.0, 0.0)

# Define a point
point = ot.Point(distribution.getDimension(), 5.0)
ott.assert_almost_equal(point, [5.0], 0.0, 0.0)

# PDF value
PDF = distribution.computePDF(point)
ott.assert_almost_equal(PDF, 0.002980, 1e-6, 1e-6)

# derivative of the PDF with regards the parameters of the distribution
CDF = distribution.computeCDF(point)
ott.assert_almost_equal(CDF, 0.003653, 1e-6, 1e-6)
# quantile
quantile = distribution.computeQuantile(0.95)
ott.assert_almost_equal(quantile, [13.0], 0.0, 0.0)
ott.assert_almost_equal(distribution.computeCDF(quantile), 0.964732, 1e-6, 1e-6)
ott.assert_almost_equal(distribution.computeEntropy(), 1.986819, 1e-6, 1e-6)
mean = distribution.getMean()
ott.assert_almost_equal(mean, [10.5], 1e-5, 0.0)
standardDeviation = distribution.getStandardDeviation()
ott.assert_almost_equal(standardDeviation, [1.77482], 1e-5, 0.0)
skewness = distribution.getSkewness()
ott.assert_almost_equal(skewness, [-0.225374], 1e-5, 0.0)
kurtosis = distribution.getKurtosis()
ott.assert_almost_equal(kurtosis, [2.91746], 1e-5, 0.0)
covariance = distribution.getCovariance()
ott.assert_almost_equal(covariance, ot.CovarianceMatrix(1, [3.15]), 1e-5, 0.0)
parameters = distribution.getParametersCollection()
assert len(parameters) == 1
ott.assert_almost_equal(parameters[0], [15.0, 0.7], 0.0, 0.0)
assert distribution.getStandardRepresentative().getImplementation().getClassName() == "Binomial"
# Confidence interval
alpha = 0.05
bounds = distribution.computeBilateralConfidenceInterval(1 - alpha)
ott.assert_almost_equal(bounds.getLowerBound(), [7.0], 0.0, 0.0)
ott.assert_almost_equal(bounds.getUpperBound(), [13.0], 0.0, 0.0)

# check survival at upper bound
distribution = ot.Binomial(10, 1.0)
assert distribution.computeCDF(10.0) == 1.0
assert distribution.computeComplementaryCDF(10.0) == 0.0
assert distribution.computeSurvivalFunction(10.0) == 0.0

# negative quantile bug
distribution = ot.Binomial(3, 0.5)
assert distribution.computeScalarQuantile(0.9, True) == 0

# degenerate distributions return their single atom at every level, see #2605
assert ot.Binomial(3, 0.0).computeScalarQuantile(0.0, True) == 0.0
assert ot.Binomial(3, 0.0).computeScalarQuantile(0.0, False) == 0.0
assert ot.Binomial(3, 0.0).computeScalarQuantile(0.7, True) == 0.0
assert ot.Binomial(3, 1.0).computeScalarQuantile(0.0, False) == 3.0
assert ot.Binomial(3, 1.0).computeScalarQuantile(0.0, True) == 3.0
assert ot.Binomial(3, 1.0).computeScalarQuantile(0.7, False) == 3.0

# first parameter range excludes 2^64 and NaN
degenerate = ot.Binomial(3, 0.5)
with ott.assert_raises(TypeError):
    degenerate.setParameter([1.8446744073709552e19, 0.5])
with ott.assert_raises(TypeError):
    degenerate.setParameter([float("nan"), 0.5])

separator = ","
path = os.path.join(os.path.dirname(__file__), "t_Binomial_std.csv")
sample = ot.Sample.ImportFromTextFile(path, separator)

rtol = 1.0e-12
atol = ot.SpecFunc.MinScalar
sample_size = sample.getSize()
for i in range(sample_size):
    x, n, pr, expected_pdf, expected_cdfp, expected_cdfq = sample[i]
    x = int(x)
    n = int(n)
    distribution = ot.Binomial(n, pr)
    computed_p = distribution.computePDF(x)
    computed_cdf = distribution.computeCDF(x)
    computed_ccdf = distribution.computeComplementaryCDF(x)
    computed_surv = distribution.computeSurvivalFunction(x)
    # Check PDF
    ott.assert_almost_equal(computed_p, expected_pdf, rtol, atol)
    # Check CDF
    ott.assert_almost_equal(computed_cdf, expected_cdfp, rtol, atol)
    # Check complementary CDF
    ott.assert_almost_equal(computed_ccdf, expected_cdfq, rtol, atol)
    # Check Survival
    ott.assert_almost_equal(computed_surv, computed_ccdf, rtol, atol)
    # Check quantile
    if pr == 0.0 or pr == 1.0:
        # Degenerate distribution: not invertible, the quantile is the
        # single atom everywhere (see dedicated asserts above)
        continue
    elif expected_cdfp < expected_cdfq:
        computed_x = int(distribution.computeQuantile(computed_cdf)[0])
        cdfXM1 = distribution.computeCDF(computed_x - 1)
        cdfX = distribution.computeCDF(computed_x)
        if computed_cdf == 0.0:
            assert computed_x == 0.0
        else:
            assert cdfXM1 < computed_cdf and computed_cdf <= cdfX
    else:
        computed_x = int(distribution.computeQuantile(computed_ccdf, True)[0])
        assert x == computed_x
    if i > 0:
        ott.assert_almost_equal(
            distribution.computeScalarQuantile(expected_cdfp), x, 0.0, 1.0
        )  # Can be off by 1 unit


# quantile bug
alpha = 0.05
beta = 0.05
for n in range(59, 100):
    d = ot.Binomial(n, alpha)
    k = d.computeQuantile(beta)[0]
    p1 = d.computeCDF(k - 1)
    p2 = d.computeCDF(k)
    ok = p1 < beta <= p2
    assert ok

ot.Log.Show(ot.Log.TRACE)
validation = ott.DistributionValidation(distribution)
validation.run()

# Corner cases

# 2147483647 is the maximum integer.
# Values greater than this are not doubles anymore.
binomial = ot.Binomial()
N = 2147483647
binomial.setN(N)
binomial.setP(1.0 / N)
computed = binomial.computePDF(1.0)
expected = 0.3678794
ott.assert_almost_equal(computed, expected, 1.0e-6, 0.0)

computed = binomial.computePDF(2.0)
expected = 0.1839397
ott.assert_almost_equal(computed, expected, 1.0e-6, 0.0)

# Extreme inputs
binomial = ot.Binomial()
binomial.setN(9999)
binomial.setP(0.5)
computed = binomial.computePDF(4999.0)
expected = 0.0079786461393821558191
ott.assert_almost_equal(computed, expected, 1.0e-7, 0.0)

# Check pdf for values of P closer to 1
binomial = ot.Binomial()
binomial.setN(2)
binomial.setP(0.9999)
computed = binomial.computePDF(1.0)
expected = 1.999799999999779835e-04
ott.assert_almost_equal(computed, expected, 1e-12, 0.0)
