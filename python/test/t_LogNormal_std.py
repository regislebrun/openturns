#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott

ot.TESTPREAMBLE()


# Instantiate one distribution object
distribution = ot.LogNormal(-1.0, 1.5, -0.5)
assert distribution.getClassName() == "LogNormal"

# Is this distribution elliptical ?
assert not distribution.isElliptical()

# Is this distribution continuous ?
assert distribution.isContinuous()

# Test for realization of distribution
oneRealization = distribution.getRealization()
ott.assert_almost_equal(oneRealization, [0.416038], 1e-5, 0.0)

# Define a point
point = ot.Point(distribution.getDimension(), 1.0)
ott.assert_almost_equal(point, [1.0], 0.0, 0.0)

# Show PDF and CDF at point
eps = 1e-5
# derivative of PDF with regards its arguments
DDF = distribution.computeDDF(point)
ott.assert_almost_equal(DDF, [-0.12381], 1e-5, 0.0)

# PDF value
LPDF = distribution.computeLogPDF(point)
ott.assert_almost_equal(LPDF, -2.168831, 1e-6, 1e-6)
PDF = distribution.computePDF(point)
ott.assert_almost_equal(PDF, 0.114311, 1e-6, 1e-6)

# derivative of the PDF with regards the parameters of the distribution
CDF = distribution.computeCDF(point)
ott.assert_almost_equal(CDF, 0.825615, 1e-6, 1e-6)
CCDF = distribution.computeComplementaryCDF(point)
ott.assert_almost_equal(CCDF, 0.174385, 1e-6, 1e-6)
CF = distribution.computeCharacteristicFunction(point[0])
ott.assert_almost_equal([CF.real, CF.imag], [0.74087680493, -0.0157786395661], 1e-9, 0.0)
PDFgr = distribution.computePDFGradient(point)
ott.assert_almost_equal(PDFgr, [0.0714046, -0.00930299, 0.12381], 1e-5, 0.0)

# derivative of the PDF with regards the parameters of the distribution
CDFgr = distribution.computeCDFGradient(point)
ott.assert_almost_equal(CDFgr, [-0.171467, -0.16066, -0.114311], 1e-5, 0.0)

# quantile
quantile = distribution.computeQuantile(0.95)
ott.assert_almost_equal(quantile, [3.83742], 1e-5, 0.0)
ott.assert_almost_equal(distribution.computeCDF(quantile), 0.95, 1e-6, 1e-6)
# Get 95% survival function
inverseSurvival = ot.Point(distribution.computeInverseSurvivalFunction(0.95))
ott.assert_almost_equal(inverseSurvival, [-0.468798], 1e-5, 0.0)
ott.assert_almost_equal(
    distribution.computeSurvivalFunction(inverseSurvival), 0.95, 1e-6, 1e-6
)
ott.assert_almost_equal(distribution.computeEntropy(), 0.824404, 1e-6, 1e-6)

# Confidence regions
interval, threshold = distribution.computeMinimumVolumeIntervalWithMarginalProbability(
    0.95
)
ott.assert_almost_equal(interval.getLowerBound(), [-0.499653], 1e-5, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [3.83753], 1e-5, 0.0)
ott.assert_almost_equal(threshold, 0.95, 1e-5, 0.0)
levelSet, beta = distribution.computeMinimumVolumeLevelSetWithThreshold(0.95)
ott.assert_almost_equal(beta, [0.0158512], 1e-5, 0.0)
interval, beta = distribution.computeBilateralConfidenceIntervalWithMarginalProbability(
    0.95
)
ott.assert_almost_equal(interval.getLowerBound(), [-0.480551], 1e-5, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [6.45838], 1e-5, 0.0)
ott.assert_almost_equal(beta, 0.95, 1e-5, 0.0)
(
    interval,
    beta,
) = distribution.computeUnilateralConfidenceIntervalWithMarginalProbability(0.95, False)
ott.assert_almost_equal(interval.getLowerBound(), [-0.5], 1e-5, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [3.83742], 1e-5, 0.0)
ott.assert_almost_equal(beta, 0.95, 1e-5, 0.0)
(
    interval,
    beta,
) = distribution.computeUnilateralConfidenceIntervalWithMarginalProbability(0.95, True)
ott.assert_almost_equal(interval.getLowerBound(), [-0.468798], 1e-5, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [35451.8], 1e-5, 0.0)
ott.assert_almost_equal(beta, 0.95, 1e-5, 0.0)

mean = distribution.getMean()
ott.assert_almost_equal(mean, [0.633148], 1e-5, 0.0)
standardDeviation = distribution.getStandardDeviation()
ott.assert_almost_equal(standardDeviation, [3.30128], 1e-5, 0.0)
skewness = distribution.getSkewness()
ott.assert_almost_equal(skewness, [33.468], 1e-5, 0.0)
kurtosis = distribution.getKurtosis()
ott.assert_almost_equal(kurtosis, [10078.3], 1e-5, 0.0)
covariance = distribution.getCovariance()
ott.assert_almost_equal(covariance, ot.CovarianceMatrix(1, [10.8985]), 1e-5, 0.0)
parameters = distribution.getParametersCollection()
assert len(parameters) == 1
ott.assert_almost_equal(parameters[0], [-1.0, 1.5, -0.5], 0.0, 0.0)
assert distribution.getStandardRepresentative().getImplementation().getClassName() == "LogNormal"

ot.Log.Show(ot.Log.TRACE)
validation = ott.DistributionValidation(distribution)
validation.setStandardDeviationTolerance(0.1)  # estimator is slow to converge
validation.setSkewnessTolerance(0.45)  # estimator is slow to converge
validation.setCorrelationTolerance(0.1)
validation.run()

# Test that setSigmaLog updates the normalization factor
dist = ot.LogNormal(0.0, 1.0, 0.0)
point = [1.0]
dist.setSigmaLog(2.0)
dist_fresh = ot.LogNormal(0.0, 2.0, 0.0)
ott.assert_almost_equal(dist.computePDF(point), dist_fresh.computePDF(point))

# NaN location parameters are rejected in all native parameter sets
with ott.assert_raises(TypeError):
    ot.LogNormalMuSigma(float("nan"), 3.3, -0.5)
with ott.assert_raises(TypeError):
    ot.LogNormalMuSigmaOverMu(float("nan"), 5.24, -0.5)
with ott.assert_raises(TypeError):
    ot.LogNormalMuErrorFactor(float("nan"), 1.5, -0.5)
