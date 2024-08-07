#! /usr/bin/env python

import openturns as ot

ot.TESTPREAMBLE()

# Instantiate one distribution object
distribution = ot.Binomial(15, 0.7)
print("Distribution ", repr(distribution))
print("Distribution ", distribution)

# Is this distribution elliptical ?
print("Elliptical = ", distribution.isElliptical())

# Is this distribution continuous ?
print("Continuous = ", distribution.isContinuous())

# Test for realization of distribution
oneRealization = distribution.getRealization()
print("oneRealization=", repr(oneRealization))

# Test for sampling
size = 10000
oneSample = distribution.getSample(size)
print("oneSample first=", repr(oneSample[0]), " last=", repr(oneSample[1]))
print("mean=", repr(oneSample.computeMean()))
print("covariance=", repr(oneSample.computeCovariance()))

# Define a point
point = ot.Point(distribution.getDimension(), 5.0)
print("Point= ", repr(point))

# Show PDF and CDF of point
eps = 1e-5
# PDF value
PDF = distribution.computePDF(point)
print("pdf     =%.6f" % PDF)
# by the finite difference technique from CDF
print(
    "pdf (FD)=%.6f"
    % (
        distribution.computeCDF(point + ot.Point(1, 0))
        - distribution.computeCDF(point + ot.Point(1, -1))
    )
)

# derivative of the PDF with regards the parameters of the distribution
CDF = distribution.computeCDF(point)
print("cdf=%.6f" % CDF)
# quantile
quantile = distribution.computeQuantile(0.95)
print("quantile=", repr(quantile))
print("cdf(quantile)=%.6f" % distribution.computeCDF(quantile))
print("entropy=%.6f" % distribution.computeEntropy())
mean = distribution.getMean()
print("mean=", repr(mean))
standardDeviation = distribution.getStandardDeviation()
print("standard deviation=", repr(standardDeviation))
skewness = distribution.getSkewness()
print("skewness=", repr(skewness))
kurtosis = distribution.getKurtosis()
print("kurtosis=", repr(kurtosis))
covariance = distribution.getCovariance()
print("covariance=", repr(covariance))
parameters = distribution.getParametersCollection()
print("parameters=", repr(parameters))
print("Standard representative=", distribution.getStandardRepresentative())
# Confidence interval
alpha = 0.05
bounds = distribution.computeBilateralConfidenceInterval(1 - alpha)
print("%.2f%% bilateral confidence interval" % ((1 - alpha) * 100), " =", bounds)

# check survival at upper bound
distribution = ot.Binomial(10, 1.0)
assert distribution.computeCDF(10.0) == 1.0
assert distribution.computeComplementaryCDF(10.0) == 0.0
assert distribution.computeSurvivalFunction(10.0) == 0.0

# negative quantile bug
distribution = ot.Binomial(3, 0.5)
assert distribution.computeScalarQuantile(0.9, True) == 0

# quantile bug
alpha = 0.05
beta = 0.05
for n in range(59, 100):
    d = ot.Binomial(n, alpha)
    k = d.computeQuantile(beta)[0]
    p1 = d.computeCDF(k - 1)
    p2 = d.computeCDF(k)
    ok = p1 < beta <= p2
    print(f"n={n} k={k:.0f} p(k-1)={p1:.4f} p(k)={p2:.4f} ok={ok}")
    assert ok
