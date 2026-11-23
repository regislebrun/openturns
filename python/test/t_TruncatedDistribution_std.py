#! /usr/bin/env python

import math as m

import openturns as ot
import openturns.experimental as otexp
import openturns.testing as ott

# Truncation of a univariate normal: the truncated distribution must agree with
# the closed-form TruncatedNormal reference
univariateCases = [
    (
        ot.TruncatedDistribution(ot.Normal(2.0, 1.5), ot.Interval(1.0, 4.0)),
        ot.TruncatedNormal(2.0, 1.5, 1.0, 4.0),
    ),
    (
        ot.TruncatedDistribution(
            ot.Normal(2.0, 1.5), 1.0, ot.TruncatedDistribution.LOWER
        ),
        ot.TruncatedNormal(2.0, 1.5, 1.0, 200.0),
    ),
    (
        ot.TruncatedDistribution(
            ot.Normal(2.0, 1.5), 4.0, ot.TruncatedDistribution.UPPER
        ),
        ot.TruncatedNormal(2.0, 1.5, -200.0, 4.0),
    ),
    (
        ot.TruncatedDistribution(
            ot.Normal(2.0, 1.5), ot.Interval([1.0], [4.0], [True], [True])
        ),
        ot.TruncatedNormal(2.0, 1.5, 1.0, 4.0),
    ),
]

# Truncation of a multivariate normal: there is no closed-form reference, the
# moments are checked at the end of this file
normal2d = ot.Normal(2)
bounds2d = ot.Interval([-0.5] * 2, [2.0] * 2)
truncated2d = ot.TruncatedDistribution(normal2d, bounds2d)

# Truncating a distribution whose support is already inside the interval must
# leave it unchanged
weibull = ot.WeibullMin(2.0, 3.0)

ot.RandomGenerator.SetSeed(0)
for currentDistribution, referenceDistribution in (
    univariateCases + [(truncated2d, None), (weibull, weibull)]
):
    dimension = currentDistribution.getDimension()
    hasReference = referenceDistribution is not None
    # a point inside every truncation interval
    point = ot.Point(dimension, 1.5)

    # A truncated normal is not elliptical, and every case is continuous
    assert not currentDistribution.isElliptical()
    assert currentDistribution.isContinuous()

    # The realizations stay in the truncation interval, and the sample mean
    # converges to the theoretical one
    oneRealization = currentDistribution.getRealization()
    assert currentDistribution.getRange().contains(oneRealization)
    sampleSize = 10000
    sampleMean = currentDistribution.getSample(sampleSize).computeMean()
    mean = currentDistribution.getMean()
    standardDeviation = currentDistribution.getStandardDeviation()
    for i in range(dimension):
        assert abs(sampleMean[i] - mean[i]) < 5.0 * standardDeviation[i] / m.sqrt(
            sampleSize
        )

    # Density, first derivative and CDF against the closed-form reference, or
    # against the normalization by the probability of the truncation box
    if hasReference:
        ott.assert_almost_equal(
            currentDistribution.computeDDF(point),
            referenceDistribution.computeDDF(point),
            1e-10,
            1e-12,
        )
        ott.assert_almost_equal(
            currentDistribution.computePDF(point),
            referenceDistribution.computePDF(point),
            1e-10,
            1e-12,
        )
        ott.assert_almost_equal(
            currentDistribution.computeCDF(point),
            referenceDistribution.computeCDF(point),
            1e-10,
            1e-12,
        )
    else:
        ott.assert_almost_equal(
            currentDistribution.computePDF(point)
            * normal2d.computeProbability(bounds2d),
            normal2d.computePDF(point),
            1e-10,
            1e-12,
        )
    ott.assert_almost_equal(
        currentDistribution.computeCDF(point)
        + currentDistribution.computeComplementaryCDF(point),
        1.0,
        1e-10,
        1e-12,
    )

    # The parameter vector and its description must have the same size: an
    # infinite bound is not part of the parameter vector
    parameter = currentDistribution.getParameter()
    assert (
        currentDistribution.getParameterDescription().getSize()
        == parameter.getDimension()
    )
    assert len(currentDistribution.getParametersCollection()) == 1

    # The marginal distribution of a truncated distribution is the marginal of
    # the truncated distribution, not the truncation of the marginal
    marginal = currentDistribution.getMarginal(0)
    if dimension == 1:
        ott.assert_almost_equal(
            marginal.computePDF(1.5), currentDistribution.computePDF(1.5), 1e-12, 1e-14
        )
    else:
        assert marginal.getImplementation().getClassName() == "TruncatedDistribution"
        assert marginal.computeCDF([1.5]) > 0.0
    assert currentDistribution.getStandardRepresentative().getDimension() == dimension

    # A conditional distribution without conditioning is the distribution
    # itself, and the sequential form is the same quantity
    conditionalPDF = currentDistribution.computeConditionalPDF(mean[0], [])
    sequentialPDF = currentDistribution.computeSequentialConditionalPDF(mean)
    ott.assert_almost_equal(conditionalPDF, sequentialPDF[0], 1e-12, 1e-14)
    if dimension == 1:
        ott.assert_almost_equal(
            conditionalPDF,
            currentDistribution.computePDF(ot.Point([mean[0]])),
            1e-12,
            1e-14,
        )
    conditionalCDF = currentDistribution.computeConditionalCDF(mean[0], [])
    sequentialCDF = currentDistribution.computeSequentialConditionalCDF(mean)
    ott.assert_almost_equal(conditionalCDF, sequentialCDF[0], 1e-12, 1e-14)

    # Quantile, survival function and confidence regions are only defined for
    # a univariate distribution
    if dimension == 1:
        quantile = currentDistribution.computeQuantile(0.95)
        ott.assert_almost_equal(
            currentDistribution.computeCDF(quantile), 0.95, 1e-10, 1e-12
        )
        inverseSurvival = ot.Point(
            currentDistribution.computeInverseSurvivalFunction(0.95)
        )
        ott.assert_almost_equal(
            currentDistribution.computeSurvivalFunction(inverseSurvival),
            0.95,
            1e-10,
            1e-12,
        )
        if hasReference:
            ott.assert_almost_equal(
                quantile, referenceDistribution.computeQuantile(0.95), 1e-10, 1e-12
            )
            ott.assert_almost_equal(
                inverseSurvival,
                referenceDistribution.computeInverseSurvivalFunction(0.95),
                1e-10,
                1e-12,
            )
        # Each confidence region carries the requested probability and stays
        # inside the numerical range
        numericalRange = currentDistribution.getRange()
        confidenceRegions = [
            currentDistribution.computeMinimumVolumeIntervalWithMarginalProbability(
                0.95
            ),
            currentDistribution.computeBilateralConfidenceIntervalWithMarginalProbability(
                0.95
            ),
            currentDistribution.computeUnilateralConfidenceIntervalWithMarginalProbability(
                0.95, False
            ),
            currentDistribution.computeUnilateralConfidenceIntervalWithMarginalProbability(
                0.95, True
            ),
        ]
        for region, regionProbability in confidenceRegions:
            ott.assert_almost_equal(regionProbability, 0.95, 1e-10, 1e-12)
            lowerBound = region.getLowerBound()[0]
            upperBound = region.getUpperBound()[0]
            if numericalRange.getFiniteLowerBound()[0]:
                assert lowerBound >= numericalRange.getLowerBound()[0] - 1e-12
            if numericalRange.getFiniteUpperBound()[0]:
                assert upperBound <= numericalRange.getUpperBound()[0] + 1e-12
            ott.assert_almost_equal(
                currentDistribution.computeCDF(upperBound)
                - currentDistribution.computeCDF(lowerBound),
                0.95,
                1e-6,
                1e-8,
            )
        # The entropy of a truncated normal is the entropy of its
        # simplified version
        if hasReference:
            ott.assert_almost_equal(
                currentDistribution.computeEntropy(),
                referenceDistribution.computeEntropy(),
                1e-10,
                1e-12,
            )

    # The first four moments and the covariance against the closed-form
    # reference, or the structural properties when there is none
    covariance = currentDistribution.getCovariance()
    if hasReference:
        ott.assert_almost_equal(
            mean, referenceDistribution.getMean(), 1e-10, 1e-12
        )
        ott.assert_almost_equal(
            standardDeviation,
            referenceDistribution.getStandardDeviation(),
            1e-10,
            1e-12,
        )
        ott.assert_almost_equal(
            currentDistribution.getSkewness(),
            referenceDistribution.getSkewness(),
            1e-8,
            1e-10,
        )
        ott.assert_almost_equal(
            currentDistribution.getKurtosis(),
            referenceDistribution.getKurtosis(),
            1e-8,
            1e-10,
        )
        ott.assert_almost_equal(
            covariance, referenceDistribution.getCovariance(), 1e-10, 1e-12
        )
    else:
        assert bounds2d.contains(mean)
        for i in range(dimension):
            assert covariance[i, i] > 0.0
            for j in range(i):
                assert covariance[i, j] == covariance[j, i]
                assert abs(covariance[i, j]) <= m.sqrt(
                    covariance[i, i] * covariance[j, j]
                )

    ot.Log.Show(ot.Log.TRACE)
    ot.RandomGenerator.SetSeed(1)
    validation = ott.DistributionValidation(currentDistribution)
    validation.skipMinimumVolumeLevelSet()  # slow
    validation.skipPDFAtLowerBound()  # the density is nonzero at the finite lower bound
    validation.skipTransformation()  # deterministic Rosenblatt+FD, no RNG consumed
    validation.skipConditional()  # deterministic chain, no RNG consumed
    # NOTE: do not reduce the validation sampling sizes here: the
    # hard-coded reference values below (KernelMixture mean/covariance)
    # depend on the exact RNG stream left by these validations
    validation.run()

# Check simplification
candidates = [
    ot.Normal(1.0, 2.0),
    ot.Uniform(1.0, 2.0),
    ot.Exponential(1.0, 2.0),
    ot.TruncatedDistribution(ot.WeibullMin(), ot.Interval(1.5, 7.8)),
    ot.Beta(1.5, 6.3, -1.0, 2.0),
    ot.JointDistribution([ot.Normal()] * 2),
    ot.BlockIndependentDistribution([ot.Normal(2), ot.Normal(2)]),
    ot.BlockIndependentCopula([ot.NormalCopula(2), ot.NormalCopula(2)]),
    ot.Dirichlet([0.7, 0.3]),
    ot.TruncatedDistribution(ot.Normal(0.0, 1.0), ot.Interval(-1.0, 1.0)),
]
intervals = [
    ot.Interval(-1.0, 4.0),
    ot.Interval(0.2, 2.4),
    ot.Interval(2.5, 65.0),
    ot.Interval(2.5, 6.0),
    ot.Interval(-2.5, 6.0),
    ot.Interval(2),
    ot.Interval(4),
    ot.Interval(4),
    ot.Interval(0.2, 2.4),
    ot.Interval(-1.0, 1.0),
]
simplifiedClassNames = [
    "TruncatedNormal",
    "Uniform",
    "Exponential",
    "TruncatedDistribution",
    "Beta",
    "JointDistribution",
    "BlockIndependentDistribution",
    "BlockIndependentCopula",
    "TruncatedDistribution",
    "TruncatedNormal",
]
for candidate, interval, simplifiedClassName in zip(
    candidates, intervals, simplifiedClassNames
):
    simplified = ot.TruncatedDistribution(candidate, interval).getSimplifiedVersion()
    assert simplified.getImplementation().getClassName() == simplifiedClassName

# Check that we can set the bounds independently
truncated = ot.TruncatedDistribution()
truncated.setDistribution(ot.Normal(20.0, 7.5))
truncated.setBounds(ot.Interval([0], [80], [True], [False]))
quantile = truncated.computeQuantile(0.9)
ott.assert_almost_equal(truncated.computeCDF(quantile), 0.9, 1e-10, 1e-12)
assert quantile[0] > 20.0
# Test for issue #1190
dist = ot.Normal(6.3e-19, 2.1e-19)
dist = ot.TruncatedDistribution(dist, 4.2e-19, ot.TruncatedDistribution.LOWER)

# non-finite bound bug: the lower bound of the second coordinate is infinite,
# and the CDF must stay a finite probability of the truncated law
bounds = ot.Interval([-2, -3], [2, 3.0], [True, False], [True, True])
dist = ot.TruncatedDistribution(ot.Normal(2), bounds)
normal = ot.Normal(2)
# the generic CDF integration is only accurate to a few 1e-5 here
cdfReference = (
    normal.computeCDF([2.0, -3.0]) - normal.computeCDF([-2.0, -3.0])
) / normal.computeProbability(bounds)
ott.assert_almost_equal(dist.computeCDF([3.0, -3.0]), cdfReference, 1e-4, 1e-4)

# relative range wrt quantile epsilon issue
unif = ot.Uniform(0.0, 1e12)
trunc = ot.TruncatedDistribution(unif, ot.Interval(0.25, 2.0))
quantile = trunc.computeQuantile(0.1)
ott.assert_almost_equal(trunc.computeCDF(quantile), 0.1, 1e-10, 1e-12)
assert 0.25 <= quantile[0] <= 2.0

# n-d CDF inversion
dim = 10
unif = ot.JointDistribution([ot.Uniform(0.0, 1.0)] * dim)
trunc = ot.TruncatedDistribution(unif, ot.Interval([0.0] * dim, [1e-6] * dim))
x = trunc.getRealization()

# marginal of truncated is not truncated marginal
R = ot.CorrelationMatrix(2)
R[1, 0] = 0.8
normal = ot.Normal([0.0] * 2, [1.0] * 2, R)
trunc = ot.TruncatedDistribution(normal, ot.Interval([-0.5] * 2, [1.0] * 2))
marginal0 = trunc.getMarginal(0)
ott.assert_almost_equal(marginal0.getMean(), [0.220527])

# simplification for KernelMixture in dimension >= 3
dim = 3
kernelMixture = ot.KernelMixture(
    ot.Normal(), [0.2] * dim, ot.Normal(dim).getSample(1000)
)
bounds = ot.Interval([-1.0] * dim, [0.5] * dim)
trunc = ot.TruncatedDistribution(kernelMixture, bounds)
assert trunc.getSimplifiedVersion().getImplementation().getClassName() == "Mixture"
mean = trunc.getMean()
ott.assert_almost_equal(mean, [-0.207172, -0.193878, -0.186422])
cov = trunc.getCovariance()
ott.assert_almost_equal(
    cov,
    ot.CovarianceMatrix(
        3,
        [
            0.173185,
            -0.0152049,
            0.000857313,
            -0.0152049,
            0.164721,
            -9.17347e-05,
            0.000857313,
            -9.17347e-05,
            0.170059,
        ],
    ),
)
# MultivariateUniform truncation simplification
mvu = otexp.MultivariateUniform([0.0, 1.0], [2.0, 3.0])
truncMvu = ot.TruncatedDistribution(
    mvu, ot.Interval([0.5, 1.5], [1.5, 2.5])
)
simplified = truncMvu.getSimplifiedVersion()
assert "MultivariateUniform" in str(simplified.getClassName)
ott.assert_almost_equal(simplified.getRange().getLowerBound(), [0.5, 1.5])
ott.assert_almost_equal(simplified.getRange().getUpperBound(), [1.5, 2.5])

# check simplification logic only checks numerical range
dist = ot.Gumbel(1.0, 0.0)
numericalRange = ot.Interval(
    dist.getRange().getLowerBound(), dist.getRange().getUpperBound()
)
truncated = ot.TruncatedDistribution(dist, numericalRange)
assert truncated.getSimplifiedVersion() == dist

# Bug 2: gradient methods return zero for points inside bounds_ but outside getRange()
weibull = ot.WeibullMin()
trunc_w = ot.TruncatedDistribution(weibull, ot.Interval(-10.0, 10.0))
# getRange() = [0, 10] (intersection with WeibullMin's support), bounds_ = [-10, 10]
point_out = [-5.0]  # inside bounds_, outside getRange()
ott.assert_almost_equal(trunc_w.computePDFGradient(point_out), [0.0] * 5)
ott.assert_almost_equal(trunc_w.computeCDFGradient(point_out), [0.0] * 5)

# getMarginal must preserve variable descriptions
sample = ot.Sample(
    [
        [-90.0, 42.0, 1.0, 2.0],
        [-89.0, 43.0, 1.1, 2.1],
        [-88.0, 44.0, 0.9, 1.9],
        [-87.0, 45.0, 1.2, 2.2],
    ]
)
sample.setDescription(["a1", "a3", "P1", "Q"])
prior = ot.JointDistribution(
    [
        ot.Uniform(-100.0, -80.0),
        ot.Uniform(35.0, 50.0),
        ot.Uniform(0.0, 2.0),
        ot.Uniform(1.0, 3.0),
    ]
)
prior.setDescription(sample.getDescription())
distribution = ot.KernelSmoothing().build(sample)
truncated = ot.TruncatedDistribution(distribution, prior.getRange())
marginal = truncated.getMarginal([0, 1])
assert marginal.getDescription() == ["a1", "a3"], f"Got {list(marginal.getDescription())}"

# Test computeSequentialConditional* with custom ordering
td = ot.TruncatedDistribution(ot.Normal(2), ot.Interval([-1.0] * 2, [1.0] * 2))
x = ot.Point([0.5, 0.5])
ordering = [0, 1]
# Identity ordering should match 1-arg version
ott.assert_almost_equal(td.computeSequentialConditionalCDF(x, ordering), td.computeSequentialConditionalCDF(x), 1e-5, 1e-5)
ott.assert_almost_equal(td.computeSequentialConditionalPDF(x, ordering), td.computeSequentialConditionalPDF(x), 1e-5, 1e-5)
ott.assert_almost_equal(td.computeSequentialConditionalQuantile(ot.Point([0.3, 0.7]), ordering), td.computeSequentialConditionalQuantile(ot.Point([0.3, 0.7])), 1e-5, 1e-5)
# Reversed ordering
x_rev = ot.Point([x[1], x[0]])
r_cdf = td.computeSequentialConditionalCDF(x_rev, [1, 0])
ott.assert_almost_equal(r_cdf, [td.computeSequentialConditionalCDF(x)[1], td.computeSequentialConditionalCDF(x_rev)[1]], 1e-4, 1e-4)

# PiecewiseLinearDistribution simplification, including endpoint regression
x_trunc = ot.Point([0.0, 1.0, 2.0, 3.0])
y_trunc = ot.Point([0.1, 0.5, 0.3, 0.1])
dist_trunc = otexp.PiecewiseLinearDistribution(x_trunc, y_trunc)

# Left endpoint: truncate at original lower bound
t_left = ot.TruncatedDistribution(dist_trunc, ot.Interval(0.0, 2.0))
s_left = t_left.getSimplifiedVersion()
for test_x in [0.0, 0.5, 1.0, 1.5, 2.0]:
    ott.assert_almost_equal(s_left.computePDF(test_x), t_left.computePDF(test_x))
    ott.assert_almost_equal(s_left.computeCDF(test_x), t_left.computeCDF(test_x))

# Right endpoint: truncate at original upper bound
t_right = ot.TruncatedDistribution(dist_trunc, ot.Interval(1.0, 3.0))
s_right = t_right.getSimplifiedVersion()
for test_x in [1.0, 1.5, 2.0, 2.5, 3.0]:
    ott.assert_almost_equal(s_right.computePDF(test_x), t_right.computePDF(test_x))
    ott.assert_almost_equal(s_right.computeCDF(test_x), t_right.computeCDF(test_x))

# Both endpoints: truncate at original bounds
t_both = ot.TruncatedDistribution(dist_trunc, ot.Interval(0.0, 3.0))
s_both = t_both.getSimplifiedVersion()
for test_x in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
    ott.assert_almost_equal(s_both.computePDF(test_x), t_both.computePDF(test_x))
    ott.assert_almost_equal(s_both.computeCDF(test_x), t_both.computeCDF(test_x))

# Interior truncation (not at original bounds)
t_interior = ot.TruncatedDistribution(dist_trunc, ot.Interval(0.5, 2.5))
s_interior = t_interior.getSimplifiedVersion()
for test_x in [0.5, 1.0, 1.5, 2.0, 2.5]:
    ott.assert_almost_equal(s_interior.computePDF(test_x), t_interior.computePDF(test_x))
    ott.assert_almost_equal(s_interior.computeCDF(test_x), t_interior.computeCDF(test_x))

# Truncated discrete distribution must keep the lower atom, see issue #3361
binomial = ot.Binomial(4, 0.5)
truncated_binomial = ot.TruncatedDistribution(binomial, ot.Interval(1.0, 3.0))
ott.assert_almost_equal(truncated_binomial.computePDF(1.0), 2.0 / 7.0)
ott.assert_almost_equal(truncated_binomial.computeCDF(1.0), 2.0 / 7.0)
ott.assert_almost_equal(truncated_binomial.computeCDF(1.5), 2.0 / 7.0)
ott.assert_almost_equal(truncated_binomial.computeCDF(0.9), 0.0)
ott.assert_almost_equal(truncated_binomial.computeSurvivalFunction(1.0), 5.0 / 7.0)
ott.assert_almost_equal(truncated_binomial.computeQuantile(0.25)[0], 1.0)
ott.assert_almost_equal(truncated_binomial.computeQuantile(0.0)[0], 1.0)
ott.assert_almost_equal(truncated_binomial.computeQuantile(1.0, True)[0], 1.0)

# Truncating a 1D Dirichlet must condition, not rescale, see issue #3377
truncated_dirichlet = ot.TruncatedDistribution(
    ot.Dirichlet([2.0, 3.0]), ot.Interval(0.25, 0.75)
)
ott.assert_almost_equal(truncated_dirichlet.computePDF(0.5), 24.0 / 11.0)
ott.assert_almost_equal(truncated_dirichlet.computeCDF(0.5), 109.0 / 176.0)
ott.assert_almost_equal(truncated_dirichlet.getMean()[0], 203.0 / 440.0)

# Discrete truncation retains the lower atom: the CDF gradient must match
# the derivative of the CDF, see the retained F(a-) lower bound
truncated_poisson = ot.TruncatedDistribution(ot.Poisson(2.0), ot.Interval(1.0, 3.0))
step = 1.0e-6
cdf_fd = (
    ot.TruncatedDistribution(ot.Poisson(2.0 + step), ot.Interval(1.0, 3.0)).computeCDF(1.0)
    - ot.TruncatedDistribution(ot.Poisson(2.0 - step), ot.Interval(1.0, 3.0)).computeCDF(1.0)
) / (2.0 * step)
ott.assert_almost_equal(truncated_poisson.computeCDFGradient([1.0])[0], cdf_fd, 1e-4, 0.0)
# Discrete truncation has no bound derivatives: the full PDF gradient at x=2
# is [3/128, 0, 0], without continuous boundary terms
ott.assert_almost_equal(
    truncated_poisson.computePDFGradient([2.0]), [3.0 / 128.0, 0.0, 0.0], 1e-12, 1e-12
)

# The generic characteristic function must stay bounded for a huge argument:
# the interval count used to overflow and select an intractable Gaussian
# quadrature instead of the bounded Filon rule. TruncatedDistribution
# inherits the generic implementation, unlike Uniform which is analytic.
truncated_uniform = ot.TruncatedDistribution(
    ot.Uniform(-0.5, 1.5), ot.Interval(0.0, 1.0)
)
truncated_uniform.setIntegrationNodesNumber(256)
huge_x = 2.0**63
huge_cf = truncated_uniform.computeCharacteristicFunction(huge_x)
huge_ref = (complex(m.cos(huge_x), m.sin(huge_x)) - 1.0) / complex(0.0, huge_x)
# Scale by huge_x so that a zero result is distinguishable from the reference:
# the generic Filon rule is only accurate to a few percent at this argument,
# but the scaled value must stay close to sin(huge_x) + i (1 - cos(huge_x)).
ott.assert_almost_equal(
    [huge_cf.real * huge_x, huge_cf.imag * huge_x],
    [huge_ref.real * huge_x, huge_ref.imag * huge_x],
    0.0,
    0.25,
)

# Moments of a multivariate normal distribution truncated to a box: the
# correlated case has no closed form, TruncatedDistribution estimates them with
# the weighted Genz algorithm. The reference values below are computed
# independently, with a tensor product of Gauss-Legendre rules over the box
# (300 nodes per direction in dimension 2, 60 in dimension 3, 30 in dimension 4
# and 20 in dimension 5), converged to better than 1e-10.
truncated_normal_cases = [
    {
        "mean": [0.0, 0.0],
        "sigma": [1.0, 1.0],
        "R": [[1.0, 0.4], [0.4, 1.0]],
        "lower": [-1.0, -2.0],
        "upper": [1.5, 0.8],
        "probability": 0.5939176832963081,
        "mean_reference": [0.08469492204844642, -0.26029389970717953],
        "covariance_reference": [
            [0.39942143101324556, 0.07928170324842927],
            [0.07928170324842927, 0.4357235026137426],
        ],
    },
    {
        "mean": [0.2, -0.3, 0.5],
        "sigma": [1.0, 1.5, 0.5],
        "R": [[1.0, 0.3, -0.2], [0.3, 1.0, 0.5], [-0.2, 0.5, 1.0]],
        "lower": [-1.0, -1.0, -1.5],
        "upper": [1.5, 1.0, 0.5],
        "probability": 0.18045185542546358,
        "mean_reference": [0.40226827142531923, -0.16224738523408355, 0.17662430784873628],
        "covariance_reference": [
            [0.39180986022116915, 0.04317138953644785, -0.024697967872038158],
            [0.04317138953644786, 0.29265571892123965, 0.019107217737937263],
            [-0.024697967872038158, 0.019107217737937263, 0.061379386810656175],
        ],
    },
    {
        "mean": [0.5, -0.5, 0.0, 1.0],
        "sigma": [1.0, 0.5, 2.0, 1.0],
        "R": [
            [1.0, 0.3, -0.2, 0.1],
            [0.3, 1.0, 0.5, -0.1],
            [-0.2, 0.5, 1.0, 0.2],
            [0.1, -0.1, 0.2, 1.0],
        ],
        "lower": [-1.0, -1.0, -1.0, -1.0],
        "upper": [1.5, 1.0, 1.0, 2.0],
        "probability": 0.22149132862941065,
        "mean_reference": [
            0.39186791012228145,
            -0.41837861224541556,
            0.020972986287875535,
            0.7383749794973566,
        ],
        "covariance_reference": [
            [0.3991916571658438, 0.060010476161935236, -0.022875814767693844, 0.04246341063857634],
            [0.06001047616193525, 0.11367350490774178, 0.03207823873538321, -0.043038189130826854],
            [-0.022875814767693847, 0.03207823873538321, 0.31994100242011536, 0.02434383422598973],
            [0.04246341063857634, -0.043038189130826854, 0.02434383422598973, 0.5148624666413888],
        ],
    },
    {
        "mean": [0.5, -0.5, 0.0, 1.0, -0.2],
        "sigma": [1.0, 0.5, 2.0, 1.0, 0.75],
        "R": [
            [1.0, 0.3, -0.2, 0.1, 0.2],
            [0.3, 1.0, 0.5, -0.1, 0.1],
            [-0.2, 0.5, 1.0, 0.2, -0.1],
            [0.1, -0.1, 0.2, 1.0, 0.3],
            [0.2, 0.1, -0.1, 0.3, 1.0],
        ],
        "lower": [-1.0, -1.0, -1.0, -1.0, -1.0],
        "upper": [1.5, 1.0, 1.0, 2.0, 0.5],
        "probability": 0.1557250066773261,
        "mean_reference": [
            0.3917103492763597,
            -0.4202005896333609,
            0.02101534506258823,
            0.7422896592043253,
            -0.2475956891974728,
        ],
        "covariance_reference": [
            [0.39786008719712795, 0.05811067120866324, -0.02189017962183261, 0.03537911529156529, 0.013776263007972397],
            [0.05811067120866324, 0.11217103861050448, 0.03287149229844493, -0.05014781417249198, 0.010931168982666481],
            [-0.02189017962183261, 0.03287149229844493, 0.31976842571776165, 0.02784276430124062, -0.005563626383940408],
            [0.03537911529156529, -0.05014781417249198, 0.02784276430124062, 0.5003843611098998, 0.041588127163048355],
            [0.013776263007972397, 0.010931168982666481, -0.005563626383940409, 0.041588127163048355, 0.1618100137880299],
        ],
    },
]

for case in truncated_normal_cases:
    dimension = len(case["mean"])
    correlation = ot.CorrelationMatrix(dimension)
    for i in range(dimension):
        for j in range(i):
            correlation[i, j] = case["R"][i][j]
    normal = ot.Normal(case["mean"], case["sigma"], correlation)
    bounds = ot.Interval(case["lower"], case["upper"])
    truncated = ot.TruncatedDistribution(normal, bounds)
    # The truncation interval is entirely covered by the truncated distribution
    ott.assert_almost_equal(truncated.computeCDF(bounds.getUpperBound()), 1.0, 1e-12, 1e-12)
    ott.assert_almost_equal(truncated.computeCDF(bounds.getLowerBound()), 0.0, 1e-12, 1e-12)
    # Normal::computeProbability is exact in dimensions 2 and 3, and relies on
    # the Genz algorithm otherwise
    probability_tolerance = 1.0e-4 if dimension < 4 else 1.0e-3
    ott.assert_almost_equal(
        normal.computeProbability(bounds),
        case["probability"],
        probability_tolerance,
        probability_tolerance,
    )
    # The mean and the covariance are quasi-Monte Carlo estimates
    ott.assert_almost_equal(
        truncated.getMean(), case["mean_reference"], 5.0e-3, 5.0e-3
    )
    covariance = truncated.getCovariance()
    for i in range(dimension):
        for j in range(dimension):
            ott.assert_almost_equal(
                covariance[i, j], case["covariance_reference"][i][j], 5.0e-3, 5.0e-3
            )

# The marginal distribution of a truncated distribution is the marginal of the
# truncated distribution, not the truncation of the marginal distribution: its
# mean is the corresponding coordinate of the joint mean, which itself agrees
# with the Gauss-Legendre reference, and the density is normalized by the
# probability of the truncation interval
correlation = ot.CorrelationMatrix(2)
correlation[1, 0] = 0.4
normal = ot.Normal([0.0, 0.0], [1.0, 1.0], correlation)
bounds = ot.Interval([-1.0, -2.0], [1.5, 0.8])
truncated = ot.TruncatedDistribution(normal, bounds)
mean = truncated.getMean()
ott.assert_almost_equal(truncated.getMarginal(0).getMean(), [mean[0]], 0.0, 0.0)
ott.assert_almost_equal(truncated.getMarginal(1).getMean(), [mean[1]], 0.0, 0.0)
ott.assert_almost_equal(mean, [0.08469492204844642, -0.26029389970717953], 5e-3, 5e-3)
point = [0.2, -0.5]
ott.assert_almost_equal(
    truncated.computePDF(point),
    normal.computePDF(point) / normal.computeProbability(bounds),
    1e-12,
    1e-12,
)

# Independent components: the truncation is the product of univariate
# truncated normal distributions, whose moments are analytic
normal = ot.Normal([0.0, 0.5], [1.0, 2.0], ot.CorrelationMatrix(2))
truncated = ot.TruncatedDistribution(normal, ot.Interval([-1.0, -1.0], [1.5, 2.0]))
component_0 = ot.TruncatedNormal(0.0, 1.0, -1.0, 1.5)
component_1 = ot.TruncatedNormal(0.5, 2.0, -1.0, 2.0)
ott.assert_almost_equal(
    truncated.getMean(),
    [component_0.getMean()[0], component_1.getMean()[0]],
    1e-14,
    1e-14,
)
covariance = truncated.getCovariance()
ott.assert_almost_equal(covariance[0, 0], component_0.getCovariance()[0, 0], 1e-14, 1e-14)
ott.assert_almost_equal(covariance[1, 1], component_1.getCovariance()[0, 0], 1e-14, 1e-14)
ott.assert_almost_equal(covariance[0, 1], 0.0, 0.0, 0.0)

# A zero-variance coordinate is a Dirac at its mean: it is removed before the
# correlated computation, as in Normal::computeProbability, and does not
# contribute to the covariance
correlation = ot.CorrelationMatrix(3)
correlation[2, 1] = 0.4
normal = ot.Normal([0.0, 0.0, 0.0], [0.0, 1.0, 1.0], correlation)
truncated = ot.TruncatedDistribution(
    normal, ot.Interval([-1.0, -1.0, -2.0], [1.0, 1.5, 0.8])
)
mean = truncated.getMean()
ott.assert_almost_equal(mean[0], 0.0, 0.0, 0.0)
ott.assert_almost_equal(
    [mean[1], mean[2]], [0.08469492204844642, -0.26029389970717953], 5e-3, 5e-3
)
covariance = truncated.getCovariance()
ott.assert_almost_equal(covariance[0, 0], 0.0, 0.0, 0.0)
ott.assert_almost_equal(covariance[0, 1], 0.0, 0.0, 0.0)
ott.assert_almost_equal(covariance[0, 2], 0.0, 0.0, 0.0)

# A semi-infinite truncation interval is equivalent to a truncation at ten
# standard deviations: the reference values are those of the same box with a
# lower bound equal to -10, computed with the Gauss-Legendre rule
correlation = ot.CorrelationMatrix(2)
correlation[1, 0] = 0.4
normal = ot.Normal([0.0, 0.0], [1.0, 1.0], correlation)
truncated = ot.TruncatedDistribution(
    normal, ot.Interval([-1.0, -10.0], [1.5, 0.8], [True, False], [True, True])
)
ott.assert_almost_equal(
    truncated.getMean(),
    [0.07801852113755958, -0.3006921527796557],
    5e-3,
    5e-3,
)
covariance = truncated.getCovariance()
ott.assert_almost_equal(covariance[0, 0], 0.39955656670382916, 5e-3, 5e-3)
ott.assert_almost_equal(covariance[0, 1], 0.09156947577302937, 5e-3, 5e-3)
ott.assert_almost_equal(covariance[1, 1], 0.511308547747231, 5e-3, 5e-3)

# The covariance is invariant under a translation of the mean and the box,
# while the mean is translated: the second moment is accumulated around the
# box, not around the origin
correlation = ot.CorrelationMatrix(2)
correlation[1, 0] = 0.4
normal = ot.Normal([0.0, 0.0], [1.0, 1.0], correlation)
reference = ot.TruncatedDistribution(normal, ot.Interval([-1.0, -2.0], [1.5, 0.8]))
offset = 1.0e6
shifted = ot.TruncatedDistribution(
    ot.Normal([offset, offset], [1.0, 1.0], correlation),
    ot.Interval([offset - 1.0, offset - 2.0], [offset + 1.5, offset + 0.8]),
)
shifted_mean = shifted.getMean()
ott.assert_almost_equal(
    [shifted_mean[0] - offset, shifted_mean[1] - offset],
    reference.getMean(),
    1e-8,
    1e-8,
)
shifted_covariance = shifted.getCovariance()
reference_covariance = reference.getCovariance()
for i in range(2):
    for j in range(2):
        ott.assert_almost_equal(
            shifted_covariance[i, j], reference_covariance[i, j], 1e-8, 1e-8
        )

# A tail box is estimated with the same weighted Genz pass, with a probability
# of 4.1e-5: the reference values are computed with the Gauss-Legendre rule
correlation = ot.CorrelationMatrix(2)
correlation[1, 0] = 0.4
normal = ot.Normal([0.0, 0.0], [1.0, 1.0], correlation)
truncated = ot.TruncatedDistribution(normal, ot.Interval([3.0, 3.0], [4.0, 4.0]))
ott.assert_almost_equal(
    truncated.getMean(), [3.311439261484444, 3.311439261484444], 5e-3, 5e-3
)
covariance = truncated.getCovariance()
ott.assert_almost_equal(covariance[0, 0], 0.06027468010048942, 5e-3, 5e-3)
ott.assert_almost_equal(covariance[0, 1], 0.0017452711331529258, 5e-3, 5e-3)
ott.assert_almost_equal(covariance[1, 1], 0.06027468010048764, 5e-3, 5e-3)
assert abs(covariance[0, 1]) <= m.sqrt(covariance[0, 0] * covariance[1, 1])
