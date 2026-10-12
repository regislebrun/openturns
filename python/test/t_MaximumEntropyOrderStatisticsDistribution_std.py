#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott

ot.TESTPREAMBLE()


# Force the use of the approximation to avoid timeout
ot.ResourceMap.SetAsBool(
    "MaximumEntropyOrderStatisticsDistribution-UseApproximation", True
)

# Instantiate one distribution object
distribution = ot.MaximumEntropyOrderStatisticsDistribution(
    [
        ot.Trapezoidal(-2.0, -1.1, -1.0, 1.0),
        ot.LogUniform(1.0, 1.2),
        ot.Triangular(3.0, 4.5, 5.0),
        ot.Beta(2.5, 3.5, 4.7, 5.2),
    ]
)
assert distribution.getClassName() == "MaximumEntropyOrderStatisticsDistribution"

dim = distribution.getDimension()

# Is this distribution elliptical ?
assert not distribution.isElliptical()

# Test for realization of distribution
oneRealization = distribution.getRealization()
ott.assert_almost_equal(oneRealization, [-0.514848, 3.2432, 3.67427, 4.75124], 1e-5, 0.0)

# Test for sampling
size = 10000
oneSample = distribution.getSample(size)
ott.assert_almost_equal(
    oneSample[0], [-1.01206, 3.29988, 4.7226, 4.9015], 1e-5, 0.0
)
ott.assert_almost_equal(
    oneSample[size - 1], [-0.845925, 3.08338, 3.64246, 4.96147], 1e-5, 0.0
)
ott.assert_almost_equal(
    oneSample.computeMean(), [-0.690223, 3.00869, 4.17387, 4.90804], 1e-5, 0.0
)
ott.assert_almost_equal(
    oneSample.computeCovariance(),
    ot.CovarianceMatrix(
        4,
        [
            0.403946,
            -0.00036587,
            0.000932395,
            -5.81027e-05,
            -0.00036587,
            0.0300943,
            0.00129843,
            -2.28507e-05,
            0.000932395,
            0.00129843,
            0.183603,
            0.00207535,
            -5.81027e-05,
            -2.28507e-05,
            0.00207535,
            0.00865565,
        ],
    ),
    1e-5,
    0.0,
)

# Define a point
point = ot.Point([0.0, 3.2, 4.2, 5.0])
ott.assert_almost_equal(point, [0.0, 3.2, 4.2, 5.0], 0.0, 0.0)

# derivative of PDF with regards its arguments
DDF = distribution.computeDDF(point)
ott.assert_almost_equal(DDF, [-0.645398, 0.107242, 0.537832, -4.84052], 1e-5, 0.0)
# PDF value
LPDF = distribution.computeLogPDF(point)
ott.assert_almost_equal(LPDF, -0.437888, 1e-6, 1e-6)
PDF = distribution.computePDF(point)
ott.assert_almost_equal(PDF, 0.645398, 1e-6, 1e-6)
condPDF = distribution.computeConditionalPDF(4.0, [0.0, 2.0])
ott.assert_almost_equal(condPDF, 0.0, 0.0, 1e-6)

# derivative of the PDF with regards the parameters of the distribution
CDF = distribution.computeCDF(point)
ott.assert_almost_equal(CDF, 0.272530, 1e-6, 1e-6)
CCDF = distribution.computeComplementaryCDF(point)
ott.assert_almost_equal(CCDF, 0.727470, 1e-6, 1e-6)
condCDF = distribution.computeConditionalCDF(4.0, [0.0, 2.0])
ott.assert_almost_equal(condCDF, 0.0, 0.0, 1e-6)
# quantile
quantile = distribution.computeQuantile(0.95)
ott.assert_almost_equal(quantile, [0.719851, 3.31172, 4.88749, 5.11511], 1e-5, 0.0)
ott.assert_almost_equal(distribution.computeCDF(quantile), 0.95, 1e-6, 1e-6)
condQuantile = distribution.computeConditionalQuantile(0.4, [0.0, 2.0])
ott.assert_almost_equal(condQuantile, 4.088961, 1e-6, 1e-6)

mean = distribution.getMean()
ott.assert_almost_equal(mean, [-0.689247, 3.00918, 4.16667, 4.90833], 1e-5, 0.0)
standardDeviation = distribution.getStandardDeviation()
ott.assert_almost_equal(standardDeviation, [0.627114, 0.173677, 0.424918, 0.0931695], 1e-5, 0.0)
skewness = distribution.getSkewness()
ott.assert_almost_equal(skewness, [0.33917, 0.0692634, -0.422404, 0.223607], 1e-5, 0.0)
kurtosis = distribution.getKurtosis()
ott.assert_almost_equal(kurtosis, [2.39815, 1.80525, 2.4, 2.4], 1e-5, 0.0)
ot.ResourceMap.SetAsUnsignedInteger("GaussKronrod-MaximumSubIntervals", 20)
ot.ResourceMap.SetAsScalar("GaussKronrod-MaximumError", 1.0e-4)
covariance = distribution.getCovariance()
ott.assert_almost_equal(
    covariance,
    ot.CovarianceMatrix(
        4,
        [
            0.393268,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0301637,
            0.00176556,
            0.0,
            0.0,
            0.00176556,
            0.180556,
            0.00188803,
            0.0,
            0.0,
            0.00188803,
            0.00868056,
        ],
    ),
    1e-5,
    0.0,
)
correlation = distribution.getCorrelation()
ott.assert_almost_equal(
    correlation,
    ot.CorrelationMatrix(
        4,
        [
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.023924,
            0.0,
            0.0,
            0.023924,
            1.0,
            0.0476904,
            0.0,
            0.0,
            0.0476904,
            1.0,
        ],
    ),
    1e-5,
    0.0,
)
spearman = distribution.getSpearmanCorrelation()
ott.assert_almost_equal(
    spearman,
    ot.CorrelationMatrix(
        4,
        [
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0168547,
            0.0,
            0.0,
            0.0168547,
            1.0,
            0.0488836,
            0.0,
            0.0,
            0.0488836,
            1.0,
        ],
    ),
    1e-5,
    0.0,
)
ot.ResourceMap.SetAsUnsignedInteger("GaussKronrod-MaximumSubIntervals", 100)
ot.ResourceMap.SetAsScalar("GaussKronrod-MaximumError", 1.0e-12)
parameters = distribution.getParametersCollection()
assert len(parameters) == 4
ott.assert_almost_equal(parameters[0], [-2.0, -1.1, -1.0, 1.0], 0.0, 0.0)
ott.assert_almost_equal(parameters[1], [1.0, 1.2], 0.0, 0.0)
ott.assert_almost_equal(parameters[2], [3.0, 4.5, 5.0], 0.0, 0.0)
ott.assert_almost_equal(parameters[3], [2.5, 3.5, 4.7, 5.2], 0.0, 0.0)

# Extract the marginals
expectedMargins = ["Trapezoidal", "LogUniform", "Triangular", "Beta"]
expectedMarginPDF = [0.322581, 1.562500, 0.800000, 2.554910]
expectedMarginCDF = [0.838710, 0.815754, 0.480000, 0.819685]
expectedMarginQuantile = [[-0.760682], [3.00417], [4.22474], [4.90342]]
expectedMarginRealization = [[-0.79344], [2.75513], [4.11996], [4.84722]]
expectedMarginLower = [[-2.0], [2.71828], [3.0], [4.7]]
expectedMarginUpper = [[1.0], [3.32012], [5.0], [5.2]]
for i in range(dim):
    margin = distribution.getMarginal(i)
    assert margin.getImplementation().getClassName() == expectedMargins[i]
    ott.assert_almost_equal(margin.computePDF(point[i]), expectedMarginPDF[i], 1e-6, 1e-6)
    ott.assert_almost_equal(margin.computeCDF(point[i]), expectedMarginCDF[i], 1e-6, 1e-6)
    ott.assert_almost_equal(margin.computeQuantile(0.5), expectedMarginQuantile[i], 1e-5, 0.0)
    ott.assert_almost_equal(margin.getRealization(), expectedMarginRealization[i], 1e-5, 0.0)
    ott.assert_almost_equal(margin.getRange().getLowerBound(), expectedMarginLower[i], 1e-5, 0.0)
    ott.assert_almost_equal(margin.getRange().getUpperBound(), expectedMarginUpper[i], 1e-5, 0.0)

# Extract a 2-D marginal
indices = [0, 1]
assert indices == [0, 1]
margins = distribution.getMarginal(indices)
assert margins.getImplementation().getClassName() == "MaximumEntropyOrderStatisticsDistribution"
ott.assert_almost_equal(margins.computePDF([point[i] for i in indices]), 0.504032, 1e-6, 1e-6)
ott.assert_almost_equal(margins.computeCDF([point[i] for i in indices]), 0.684181, 1e-6, 1e-6)
quantile = ot.Point(margins.computeQuantile(0.8))
ott.assert_almost_equal(quantile, [0.190956, 3.25075], 1e-5, 0.0)
ott.assert_almost_equal(margins.computeCDF(quantile), 0.8, 1e-6, 1e-6)
marginsRealization = margins.getRealization()
ott.assert_almost_equal(marginsRealization, [-1.23466, 3.19381], 1e-5, 0.0)

# with Beta(i, n-i+1) marginals the distribution is the joint order statistics of a uniform sample
collBeta = [ot.Beta(1.0, 3.0, 0.0, 1.0), ot.Beta(2.0, 2.0, 0.0, 1.0), ot.Beta(3.0, 1.0, 0.0, 1.0)]
betaDist = ot.MaximumEntropyOrderStatisticsDistribution(collBeta, False)
ott.assert_almost_equal(betaDist.computeConditionalCDF(0.5, [0.7]), 0.0)
ott.assert_almost_equal(betaDist.computeConditionalPDF(0.5, [0.2]), 1.5625)
ott.assert_almost_equal(
    betaDist.computeCDF([0.2, 0.9, 0.5]), betaDist.computeCDF([0.2, 0.5, 0.5])
)
ott.assert_almost_equal(betaDist.computeCDF([0.2, 0.5, 0.8]), 0.242)
# the marginals of the Beta case are the marginals of the uniform order statistics
unifDist = ot.UniformOrderStatistics(3)
for indices in ([0, 1], [0, 2], [1, 2]):
    betaMarginal = betaDist.getMarginal(indices)
    unifMarginal = unifDist.getMarginal(indices)
    pointBeta = [[0.2, 0.5, 0.8][i] for i in indices]
    ott.assert_almost_equal(
        betaMarginal.computeCDF(pointBeta), unifMarginal.computeCDF(pointBeta)
    )
    ott.assert_almost_equal(
        betaMarginal.computePDF(pointBeta), unifMarginal.computePDF(pointBeta)
    )

ot.Log.Show(ot.Log.TRACE)
validation = ott.DistributionValidation(distribution)
validation.skipMoments()  # slow
validation.skipCorrelation()  # slow
validation.skipParameters()  # slow
validation.skipConditional()  # FIXME
validation.skipTransformation()  # FIXME
validation.run()

# CDF grid overflow is rejected in integer arithmetic before any allocation
highDim = ot.MaximumEntropyOrderStatisticsDistribution(
    [ot.Uniform(0.0, float(i + 1)) for i in range(17)]
)
ot.ResourceMap.SetAsBool(
    "MaximumEntropyOrderStatisticsDistribution-UseApproximation", False
)
try:
    with ott.assert_raises(TypeError):
        highDim.computeCDF([0.1 * (i + 1) for i in range(17)])
finally:
    ot.ResourceMap.SetAsBool(
        "MaximumEntropyOrderStatisticsDistribution-UseApproximation", True
    )
