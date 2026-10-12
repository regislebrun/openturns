#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott

ot.TESTPREAMBLE()


# Instantiate one distribution object
dim = 2
copula = ot.AliMikhailHaqCopula(0.5)
assert copula.getClassName() == "AliMikhailHaqCopula"
assert copula.getDimension() == 2
ott.assert_almost_equal(copula.getMean(), [0.5, 0.5], 1e-5, 0.0)
ott.assert_almost_equal(
    copula.getCovariance(),
    ot.CovarianceMatrix(2, [0.0833333, 0.0160319, 0.0160319, 0.0833333]),
    1e-5,
    0.0,
)

# Is this copula an elliptical distribution?
assert not copula.isElliptical()

# Is this copula elliptical ?
assert not copula.hasEllipticalCopula()

# Is this copula independent ?
assert not copula.hasIndependentCopula()

# Test for realization of distribution
oneRealization = copula.getRealization()
ott.assert_almost_equal(oneRealization, [0.629877, 0.895761], 1e-5, 0.0)

# Test for sampling
size = 10
oneSample = copula.getSample(size)
ott.assert_almost_equal(
    oneSample,
    [
        [0.135276, 0.0211783],
        [0.347057, 0.964145],
        [0.92068, 0.606772],
        [0.0632061, 0.189429],
        [0.714382, 0.449122],
        [0.373767, 0.713037],
        [0.883503, 0.393905],
        [0.92851, 0.86981],
        [0.684575, 0.852794],
        [0.359802, 0.947852],
    ],
    1e-5,
    0.0,
)

# Define a point
point = ot.Point(dim, 0.2)

# Show PDF and CDF of point
pointPDF = copula.computePDF(point)
pointCDF = copula.computeCDF(point)
ott.assert_almost_equal(point, [0.2, 0.2], 0.0, 0.0)
ott.assert_almost_equal(pointPDF, 1.208528, 1e-6, 1e-6)
ott.assert_almost_equal(pointCDF, 0.058824, 1e-6, 1e-6)

# Get 50% quantile
quantile = copula.computeQuantile(0.5)
ott.assert_almost_equal(quantile, [0.689898, 0.689898], 1e-5, 0.0)
ott.assert_almost_equal(copula.computeCDF(quantile), 0.5, 1e-6, 1e-6)

# Get 5% quantile
quantile = copula.computeQuantile(0.95, True)
ott.assert_almost_equal(quantile, [0.182457, 0.182457], 1e-5, 0.0)

# Get 95% survival function
inverseSurvival = ot.Point(copula.computeInverseSurvivalFunction(0.95))
ott.assert_almost_equal(inverseSurvival, [0.025625, 0.025625], 1e-5, 0.0)
ott.assert_almost_equal(
    copula.computeSurvivalFunction(inverseSurvival), 0.95, 1e-6, 1e-6
)
ott.assert_almost_equal(copula.computeEntropy(), -0.019288, 1e-6, 1e-6)
# Confidence regions
interval, threshold = copula.computeMinimumVolumeIntervalWithMarginalProbability(0.95)
ott.assert_almost_equal(interval.getLowerBound(), [0.0126791, 0.0126791], 1e-5, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [0.987321, 0.987321], 1e-5, 0.0)
ott.assert_almost_equal(threshold, 0.974642, 1e-5, 0.0)
levelSet, beta = copula.computeMinimumVolumeLevelSetWithThreshold(0.95)
ott.assert_almost_equal(beta, 0.709465, 1e-5, 0.0)
interval, beta = copula.computeBilateralConfidenceIntervalWithMarginalProbability(0.95)
ott.assert_almost_equal(interval.getLowerBound(), [0.0126791, 0.0126791], 1e-5, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [0.987321, 0.987321], 1e-5, 0.0)
ott.assert_almost_equal(beta, 0.974642, 1e-5, 0.0)
interval, beta = copula.computeUnilateralConfidenceIntervalWithMarginalProbability(
    0.95, False
)
ott.assert_almost_equal(interval.getLowerBound(), [0.0, 0.0], 0.0, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [0.974521, 0.974521], 1e-5, 0.0)
ott.assert_almost_equal(beta, 0.974521, 1e-5, 0.0)
interval, beta = copula.computeUnilateralConfidenceIntervalWithMarginalProbability(
    0.95, True
)
ott.assert_almost_equal(interval.getLowerBound(), [0.025625, 0.025625], 1e-5, 0.0)
ott.assert_almost_equal(interval.getUpperBound(), [1.0, 1.0], 1e-5, 0.0)
ott.assert_almost_equal(beta, 0.974375, 1e-5, 0.0)

# Extract the marginals
for i in range(dim):
    margin = copula.getMarginal(i)
    assert margin.getImplementation().getClassName() == "IndependentCopula"
    ott.assert_almost_equal(margin.computePDF([0.25]), 1.0, 1e-6, 1e-6)
    ott.assert_almost_equal(margin.computeCDF([0.25]), 0.25, 1e-6, 1e-6)
    ott.assert_almost_equal(margin.computeQuantile(0.95), [0.95], 1e-5, 0.0)
    marginRealization = margin.getRealization()
    assert marginRealization[0] >= 0.0 and marginRealization[0] <= 1.0
    if i == 0:
        ott.assert_almost_equal(marginRealization, [0.592662], 1e-5, 0.0)
    else:
        ott.assert_almost_equal(marginRealization, [0.744515], 1e-5, 0.0)

# Extract a 2-D marginal
indices = [1, 0]
assert indices == [1, 0]
margins = copula.getMarginal(indices)
assert margins.getImplementation().getClassName() == "AliMikhailHaqCopula"
ott.assert_almost_equal(margins.computePDF([0.25] * 2), 1.136188, 1e-6, 1e-6)
ott.assert_almost_equal(margins.computeCDF([0.25] * 2), 0.086957, 1e-6, 1e-6)
quantile = ot.Point(margins.computeQuantile(0.95))
ott.assert_almost_equal(quantile, [0.974521, 0.974521], 1e-5, 0.0)
ott.assert_almost_equal(margins.computeCDF(quantile), 0.95, 1e-6, 1e-6)
marginsRealization = margins.getRealization()
ott.assert_almost_equal(marginsRealization, [0.0405473, 0.677683], 1e-5, 0.0)

ot.Log.Show(ot.Log.TRACE)
validation = ott.DistributionValidation(copula)
validation.run()

# Theta setter rejects NaN and out-of-range values, accepts endpoints
probe = ot.AliMikhailHaqCopula(0.5)
with ott.assert_raises(TypeError):
    probe.setTheta(float("nan"))
with ott.assert_raises(TypeError):
    probe.setTheta(-1.5)
with ott.assert_raises(TypeError):
    probe.setTheta(1.5)
probe.setTheta(1.0)
ott.assert_almost_equal(probe.getTheta(), 1.0, 0.0, 0.0)

# Theta=1 is the uv/(u+v-uv) copula, of generator 1/t-1
unitCopula = ot.AliMikhailHaqCopula(1.0)
ott.assert_almost_equal(unitCopula.computeCDF([0.2, 0.2]), 1.0 / 9.0, 1e-12, 0.0)
ott.assert_almost_equal(unitCopula.computePDF([0.2, 0.2]), 1.71467764, 1e-6, 0.0)
ott.assert_almost_equal(unitCopula.computeArchimedeanGenerator(0.5), 1.0, 0.0, 0.0)
ott.assert_almost_equal(unitCopula.computeInverseArchimedeanGenerator(1.0), 0.5, 0.0, 0.0)
ott.assert_almost_equal(unitCopula.computeArchimedeanGeneratorDerivative(0.5), -4.0, 0.0, 0.0)
ott.assert_almost_equal(unitCopula.computeArchimedeanGeneratorSecondDerivative(0.5), 16.0, 0.0, 0.0)
