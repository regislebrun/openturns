#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott

ot.TESTPREAMBLE()

# Instantiate one distribution object
distribution = ot.Gumbel(0.5, -0.5)
assert distribution.getClassName() == "Gumbel"

# Is this distribution elliptical ?
assert not distribution.isElliptical()

# Is this distribution continuous ?
assert distribution.isContinuous()

# Test for realization of distribution
oneRealization = distribution.getRealization()
ott.assert_almost_equal(oneRealization, [-0.114155], 1e-5, 0.0)

# Define a point
point = ot.Point(distribution.getDimension(), 1.0)
ott.assert_almost_equal(point, [1.0], 0.0, 0.0)

# Show PDF and CDF at point
eps = 1e-5
# derivative of PDF with regards its arguments
DDF = distribution.computeDDF(point)
ott.assert_almost_equal(DDF, [-0.180043], 1e-5, 0.0)

# PDF value
LPDF = distribution.computeLogPDF(point)
ott.assert_almost_equal(LPDF, -2.356640, 1e-6, 1e-6)
PDF = distribution.computePDF(point)
ott.assert_almost_equal(PDF, 0.094738, 1e-6, 1e-6)

# derivative of the PDF with regards the parameters of the distribution
CDF = distribution.computeCDF(point)
ott.assert_almost_equal(CDF, 0.951432, 1e-6, 1e-6)
CCDF = distribution.computeComplementaryCDF(point)
ott.assert_almost_equal(CCDF, 0.048568, 1e-6, 1e-6)
PDFgr = distribution.computePDFGradient(point)
ott.assert_almost_equal(PDFgr, [0.350652, 0.180043], 1e-5, 0.0)

# derivative of the PDF with regards the parameters of the distribution
CDFgr = distribution.computeCDFGradient(point)
ott.assert_almost_equal(CDFgr, [-0.284214, -0.094738], 1e-5, 0.0)

# quantile
quantile = distribution.computeQuantile(0.95)
ott.assert_almost_equal(quantile, [0.985098], 1e-5, 0.0)
ott.assert_almost_equal(distribution.computeCDF(quantile), 0.95, 1e-6, 1e-6)
# Get 95% survival function
inverseSurvival = ot.Point(distribution.computeInverseSurvivalFunction(0.95))
ott.assert_almost_equal(inverseSurvival, [-1.04859], 1e-5, 0.0)
ott.assert_almost_equal(
    distribution.computeSurvivalFunction(inverseSurvival), 0.95, 1e-6, 1e-6
)
ott.assert_almost_equal(distribution.computeEntropy(), 0.884068, 1e-6, 1e-6)

# Confidence regions
interval, threshold = distribution.computeMinimumVolumeIntervalWithMarginalProbability(
    0.95
)
ott.assert_almost_equal(interval.getLowerBound(), [-1.28067], 1e-5, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [1.08074], 1e-5, 0.0)
ott.assert_almost_equal(threshold, 0.95, 1e-5, 0.0)
levelSet, beta = distribution.computeMinimumVolumeLevelSetWithThreshold(0.95)
ott.assert_almost_equal(beta, [0.0812123], 1e-5, 0.0)
interval, beta = distribution.computeBilateralConfidenceIntervalWithMarginalProbability(
    0.95
)
ott.assert_almost_equal(interval.getLowerBound(), [-1.15266], 1e-5, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [1.33812], 1e-5, 0.0)
ott.assert_almost_equal(beta, 0.95, 1e-5, 0.0)
(
    interval,
    beta,
) = distribution.computeUnilateralConfidenceIntervalWithMarginalProbability(0.95, False)
ott.assert_almost_equal(interval.getLowerBound(), [-2.23654], 1e-5, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [0.985098], 1e-5, 0.0)
ott.assert_almost_equal(beta, 0.95, 1e-5, 0.0)
(
    interval,
    beta,
) = distribution.computeUnilateralConfidenceIntervalWithMarginalProbability(0.95, True)
ott.assert_almost_equal(interval.getLowerBound(), [-1.04859], 1e-5, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [15.6181], 1e-5, 0.0)
ott.assert_almost_equal(beta, 0.95, 1e-5, 0.0)

mean = distribution.getMean()
ott.assert_almost_equal(mean, [-0.211392], 1e-5, 0.0)
standardDeviation = distribution.getStandardDeviation()
ott.assert_almost_equal(standardDeviation, [0.641275], 1e-5, 0.0)
skewness = distribution.getSkewness()
ott.assert_almost_equal(skewness, [1.13955], 1e-5, 0.0)
kurtosis = distribution.getKurtosis()
ott.assert_almost_equal(kurtosis, [5.4], 1e-5, 0.0)
covariance = distribution.getCovariance()
ott.assert_almost_equal(covariance, ot.CovarianceMatrix(1, [0.411234]), 1e-5, 0.0)
parameters = distribution.getParametersCollection()
assert len(parameters) == 1
ott.assert_almost_equal(parameters[0], [0.5, -0.5], 0.0, 0.0)
assert distribution.getStandardRepresentative().getImplementation().getClassName() == "Gumbel"

# Test for NaN in far left tail
far_left = [-1000.0]
ott.assert_almost_equal(distribution.computePDF(far_left), 0.0)
ott.assert_almost_equal(distribution.computeDDF(far_left), [0.0])
ott.assert_almost_equal(distribution.computePDFGradient(far_left), [0.0, 0.0])
ott.assert_almost_equal(distribution.computeCDFGradient(far_left), [0.0, 0.0])

ot.Log.Show(ot.Log.TRACE)
validation = ott.DistributionValidation(distribution)
validation.run()

# Single-parameter constructor defaults the location to zero, see #2509
single = ot.Gumbel(1.0)
assert single == ot.Gumbel(1.0, 0.0)
ott.assert_almost_equal(single.getParameter(), [1.0, 0.0], 0.0, 0.0)
