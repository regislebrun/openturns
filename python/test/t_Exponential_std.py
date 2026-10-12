#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott

ot.TESTPREAMBLE()

# Instantiate one distribution object
distribution = ot.Exponential(2.5, -0.5)
assert distribution.getClassName() == "Exponential"

# Is this distribution elliptical ?
assert not distribution.isElliptical()

# Is this distribution continuous ?
assert distribution.isContinuous()

# Test for realization of distribution
oneRealization = distribution.getRealization()
ott.assert_almost_equal(oneRealization, [-0.315107], 1e-5, 0.0)

# Define a point
point = ot.Point(distribution.getDimension(), 1.0)
ott.assert_almost_equal(point, [1.0], 0.0, 0.0)

# Show PDF and CDF of point
eps = 1e-5

# derivative of PDF with regards its arguments
DDF = distribution.computeDDF(point)
ott.assert_almost_equal(DDF, [-0.146986], 1e-5, 0.0)

# PDF value
LPDF = distribution.computeLogPDF(point)
ott.assert_almost_equal(LPDF, -2.833709, 1e-6, 1e-6)
PDF = distribution.computePDF(point)
ott.assert_almost_equal(PDF, 0.058794, 1e-6, 1e-6)

# Density at the support boundary equals lambda, not zero (see #3066)
ott.assert_almost_equal(ot.Exponential().computePDF(0.0), 1.0, 0.0, 0.0)
ott.assert_almost_equal(ot.Exponential().computeLogPDF(0.0), 0.0, 0.0, 0.0)
ott.assert_almost_equal(
    ot.Exponential(2.5, -0.5).computePDF([-0.5]), 2.5, 0.0, 0.0
)

# derivative of the PDF with regards the parameters of the distribution
CDF = distribution.computeCDF(point)
ott.assert_almost_equal(CDF, 0.976482, 1e-6, 1e-6)
CCDF = distribution.computeComplementaryCDF(point)
ott.assert_almost_equal(CCDF, 0.023518, 1e-6, 1e-6)
CF = distribution.computeCharacteristicFunction(point[0])
ott.assert_almost_equal([CF.real, CF.imag], [0.921855842528, -0.110683201593], 1e-9, 0.0)
PDFgr = distribution.computePDFGradient(point)
ott.assert_almost_equal(PDFgr, [-0.0646738, 0.146986], 1e-5, 0.0)

# derivative of the PDF with regards the parameters of the distribution
CDFgr = distribution.computeCDFGradient(point)
ott.assert_almost_equal(CDFgr, [0.0352766, -0.0587944], 1e-5, 0.0)

# quantile
quantile = distribution.computeQuantile(0.95)
ott.assert_almost_equal(quantile, [0.698293], 1e-5, 0.0)
ott.assert_almost_equal(distribution.computeCDF(quantile), 0.95, 1e-6, 1e-6)
# Get 95% survival function
inverseSurvival = ot.Point(distribution.computeInverseSurvivalFunction(0.95))
ott.assert_almost_equal(inverseSurvival, [-0.479483], 1e-5, 0.0)
ott.assert_almost_equal(
    distribution.computeSurvivalFunction(inverseSurvival), 0.95, 1e-6, 1e-6
)
ott.assert_almost_equal(distribution.computeEntropy(), 0.083709, 1e-6, 1e-6)

# Confidence regions
interval, threshold = distribution.computeMinimumVolumeIntervalWithMarginalProbability(
    0.95
)
ott.assert_almost_equal(interval.getLowerBound(), [-0.5], 1e-5, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [0.698293], 1e-5, 0.0)
ott.assert_almost_equal(threshold, 0.95, 1e-5, 0.0)
levelSet, beta = distribution.computeMinimumVolumeLevelSetWithThreshold(0.95)
ott.assert_almost_equal(beta, [0.125], 1e-5, 0.0)
interval, beta = distribution.computeBilateralConfidenceIntervalWithMarginalProbability(
    0.95
)
ott.assert_almost_equal(interval.getLowerBound(), [-0.489873], 1e-5, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [0.975552], 1e-5, 0.0)
ott.assert_almost_equal(beta, 0.95, 1e-5, 0.0)
(
    interval,
    beta,
) = distribution.computeUnilateralConfidenceIntervalWithMarginalProbability(0.95, False)
ott.assert_almost_equal(interval.getLowerBound(), [-0.5], 1e-5, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [0.698293], 1e-5, 0.0)
ott.assert_almost_equal(beta, 0.95, 1e-5, 0.0)
(
    interval,
    beta,
) = distribution.computeUnilateralConfidenceIntervalWithMarginalProbability(0.95, True)
ott.assert_almost_equal(interval.getLowerBound(), [-0.479483], 1e-5, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [12.3945], 1e-5, 0.0)
ott.assert_almost_equal(beta, 0.95, 1e-5, 0.0)

mean = distribution.getMean()
ott.assert_almost_equal(mean, [-0.1], 1e-5, 0.0)
standardDeviation = distribution.getStandardDeviation()
ott.assert_almost_equal(standardDeviation, [0.4], 1e-5, 0.0)
skewness = distribution.getSkewness()
ott.assert_almost_equal(skewness, [2.0], 1e-5, 0.0)
kurtosis = distribution.getKurtosis()
ott.assert_almost_equal(kurtosis, [9.0], 1e-5, 0.0)
covariance = distribution.getCovariance()
ott.assert_almost_equal(covariance, ot.CovarianceMatrix(1, [0.16]), 1e-5, 0.0)
parameters = distribution.getParametersCollection()
assert len(parameters) == 1
ott.assert_almost_equal(parameters[0], [2.5, -0.5], 0.0, 0.0)
assert distribution.getStandardRepresentative().getImplementation().getClassName() == "Gamma"

ot.Log.Show(ot.Log.TRACE)
validation = ott.DistributionValidation(distribution)
validation.skipPDFAtLowerBound()  # the density is lambda at the finite lower bound, see #3066
validation.run()

# Right-sided derivatives at the support boundary
edge = ot.Exponential(2.5, -0.5)
ott.assert_almost_equal(edge.computePDFGradient([-0.5]), [1.0, 6.25], 1e-12, 0.0)
ott.assert_almost_equal(edge.computeDDF([-0.5]), [-6.25], 1e-12, 0.0)
