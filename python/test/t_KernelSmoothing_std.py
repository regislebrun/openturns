#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import os

ot.TESTPREAMBLE()


# Instantiate one distribution object
dim = 2
meanPoint = [0.5, -0.5]
sigma = [2.0, 3.0]
R = ot.CorrelationMatrix(dim)
for i in range(1, dim):
    R[i, i - 1] = 0.5

distribution = ot.Normal(meanPoint, sigma, R)
discretization = 100
kernel = ot.KernelSmoothing()
sample = distribution.getSample(discretization)
kernels = ot.DistributionCollection(0)
kernels.add(ot.Normal())
kernels.add(ot.Epanechnikov())
kernels.add(ot.Uniform())
kernels.add(ot.Triangular())
kernels.add(ot.Logistic())
kernels.add(ot.Beta(2.0, 2.0, -1.0, 1.0))
kernels.add(ot.Beta(3.0, 3.0, -1.0, 1.0))
expectedBandwidth = [
    [0.745017, 1.51724],
    [1.66591, 3.39265],
    [1.29041, 2.62794],
    [1.82491, 3.71646],
    [0.41075, 0.836498],
    [1.66591, 3.39265],
    [1.97113, 4.01424],
]
expectedSmoothedPDF = [0.0247623, 0.0245905, 0.0235911, 0.0246934, 0.0254113, 0.0245905, 0.0244981]
expectedSmoothedCDF = [0.279515, 0.278675, 0.278627, 0.278851, 0.280123, 0.278675, 0.27883]
expectedMixturePDF = [0.0247621, 0.0245902, 0.0235911, 0.0246934, 0.0254111, 0.0245902, 0.024498]
expectedMixtureCDF = [0.279514, 0.278675, 0.278627, 0.27885, 0.280123, 0.278675, 0.27883]
for ki, kernel in enumerate(kernels):
    smoother = ot.KernelSmoothing(kernel)
    smoothed = smoother.build(sample)
    bw = smoother.getBandwidth()
    ott.assert_almost_equal(bw, expectedBandwidth[ki], 1e-5, 0.0)
    meanSmoothed = smoothed.getMean()
    ott.assert_almost_equal(meanSmoothed, [0.50406, -0.708164], 1e-5, 0.0)
    covarianceSmoothed = smoothed.getCovariance()
    ott.assert_almost_equal(
        covarianceSmoothed,
        ot.CovarianceMatrix(2, [4.10693, 2.78725, 2.78725, 11.7692]),
        1e-5,
        0.0,
    )
    # Define a point
    point = ot.Point(smoothed.getDimension(), 0.0)

    # Show PDF and CDF of point point
    pointPDF = smoothed.computePDF(point)
    pointCDF = smoothed.computeCDF(point)
    ott.assert_almost_equal(point, [0.0, 0.0], 0.0, 0.0)
    ott.assert_almost_equal(pointPDF, expectedSmoothedPDF[ki], 1e-5, 0.0)
    ott.assert_almost_equal(distribution.computePDF(point), 0.0280502, 1e-5, 0.0)
    ott.assert_almost_equal(pointCDF, expectedSmoothedCDF[ki], 1e-5, 0.0)
    ott.assert_almost_equal(distribution.computeCDF(point), 0.305529, 1e-5, 0.0)

    # as mixture
    mixture = smoother.buildAsMixture(
        sample, smoother.computeSilvermanBandwidth(sample)
    )
    pointPDF = mixture.computePDF(point)
    pointCDF = mixture.computeCDF(point)
    ott.assert_almost_equal(pointPDF, expectedMixturePDF[ki], 1e-5, 0.0)
    ott.assert_almost_equal(distribution.computePDF(point), 0.0280502, 1e-5, 0.0)
    ott.assert_almost_equal(pointCDF, expectedMixtureCDF[ki], 1e-5, 0.0)
    ott.assert_almost_equal(distribution.computeCDF(point), 0.305529, 1e-5, 0.0)


# Test for boundary correction
distributionCollection = ot.DistributionCollection(2)
distributionCollection[0] = ot.Normal(0.0, 1.0)
distributionCollection[1] = ot.Beta(0.7, 0.9, -1.0, 2.0)
sampleCollection = [
    distributionCollection[0].getSample(discretization),
    distributionCollection[1].getSample(discretization),
]
expectedBoundaryPDF = [
    [[0.302612, 0.302613], [0.339484, 0.537733]],
    [[0.302157, 0.302157], [0.316357, 0.521347]],
    [[0.300471, 0.300471], [0.292855, 0.491997]],
    [[0.302471, 0.302471], [0.334081, 0.531359]],
    [[0.303756, 0.303774], [0.355609, 0.546998]],
    [[0.302157, 0.302157], [0.316357, 0.521347]],
    [[0.302041, 0.302041], [0.324187, 0.527483]],
]
expectedBoundaryCDF = [
    [[0.206809, 0.206809], [0.0881112, 0.0535112]],
    [[0.208576, 0.208576], [0.0891315, 0.0516683]],
    [[0.209151, 0.209151], [0.0900747, 0.0507867]],
    [[0.208301, 0.208301], [0.0887193, 0.0527025]],
    [[0.205753, 0.20575], [0.0871653, 0.05477]],
    [[0.208576, 0.208576], [0.0891315, 0.0516683]],
    [[0.207977, 0.207977], [0.0888216, 0.0522068]],
]
expectedExactPDF = [[0.266085, 0.597681], [0.266085, 0.597681], [0.266085, 0.597681], [0.266085, 0.597681], [0.266085, 0.597681], [0.266085, 0.597681], [0.266085, 0.597681]]
expectedExactCDF = [[0.18406, 0.0852122], [0.18406, 0.0852122], [0.18406, 0.0852122], [0.18406, 0.0852122], [0.18406, 0.0852122], [0.18406, 0.0852122], [0.18406, 0.0852122]]
for ki in range(kernels.getSize()):
    kernel = kernels[ki]
    smoother = ot.KernelSmoothing(kernel)
    for j in range(2):
        for ci, corr in enumerate([False, True]):
            smoother.setBoundingOption(
                ot.KernelSmoothing.BOTH if corr else ot.KernelSmoothing.NONE
            )
            assert smoother.getBoundingOption() == (
                ot.KernelSmoothing.BOTH if corr else ot.KernelSmoothing.NONE
            ), "wrong bounding option"
            smoothed = smoother.build(sampleCollection[j])
            # Define a point
            point = ot.Point(smoothed.getDimension(), -0.9)

            # Show PDF and CDF of point point
            pointPDF = smoothed.computePDF(point)
            pointCDF = smoothed.computeCDF(point)
            ott.assert_almost_equal(pointPDF, expectedBoundaryPDF[ki][j][ci], 1e-5, 0.0)
            ott.assert_almost_equal(
                distributionCollection[j].computePDF(point), expectedExactPDF[ki][j], 1e-5, 0.0
            )
            ott.assert_almost_equal(pointCDF, expectedBoundaryCDF[ki][j][ci], 1e-5, 0.0)
            ott.assert_almost_equal(
                distributionCollection[j].computeCDF(point), expectedExactCDF[ki][j], 1e-5, 0.0
            )

# bounding option accessor round-trips all values, default is NONE
smoother = ot.KernelSmoothing()
assert smoother.getBoundingOption() == ot.KernelSmoothing.NONE, "wrong default"
for option in [
    ot.KernelSmoothing.NONE,
    ot.KernelSmoothing.LOWER,
    ot.KernelSmoothing.UPPER,
    ot.KernelSmoothing.BOTH,
]:
    smoother.setBoundingOption(option)
    assert smoother.getBoundingOption() == option, "wrong bounding option"

sample = ot.Normal().getSample(5000)
ks1 = ot.KernelSmoothing(ot.Normal(), True, 64).build(sample)
ks2 = ot.KernelSmoothing(ot.Normal(), True, 1024).build(sample)
ks3 = ot.KernelSmoothing(ot.Normal(), False).build(sample)
point = 0.3
ott.assert_almost_equal(ks1.computePDF(point), 0.365181, 1e-5, 0.0)
ott.assert_almost_equal(ks2.computePDF(point), 0.365423, 1e-5, 0.0)
ott.assert_almost_equal(ks3.computePDF(point), 0.365424, 1e-5, 0.0)
# Test with varying boundary corrections
coll = [ot.Uniform(-1.0, 1.0), ot.Uniform(0.0, 2.0)]
left = [-0.9, 0.1]
right = [0.9, 1.9]
expectedLeft = [
    [0.38087, 0.521626, 0.506571, 0.38087, 0.38087, 0.521626, 0.506571, 0.521626, 0.506571],
    [0.379383, 0.485563, 0.466624, 0.379383, 0.379383, 0.485563, 0.466624, 0.485563, 0.466624],
]
expectedRight = [
    [0.361147, 0.361147, 0.361147, 0.506744, 0.505152, 0.506744, 0.506744, 0.505152, 0.505152],
    [0.404847, 0.404847, 0.404847, 0.512248, 0.511777, 0.512248, 0.512248, 0.511777, 0.511777],
]
for nDist in range(len(coll)):
    baseDist = coll[nDist]
    sample = baseDist.getSample(500)
    algo1 = ot.KernelSmoothing(ot.Normal(), False)
    algo1.setBoundingOption(ot.KernelSmoothing.NONE)
    ks1 = algo1.build(sample)
    ott.assert_almost_equal(ks1.computePDF(left[nDist]), expectedLeft[nDist][0], 1e-5, 0.0)
    ott.assert_almost_equal(ks1.computePDF(right[nDist]), expectedRight[nDist][0], 1e-5, 0.0)

    algo2 = ot.KernelSmoothing(ot.Normal(), False)
    algo2.setBoundingOption(ot.KernelSmoothing.LOWER)
    algo2.setAutomaticLowerBound(True)
    ks2 = algo2.build(sample)
    ott.assert_almost_equal(ks2.computePDF(left[nDist]), expectedLeft[nDist][1], 1e-5, 0.0)
    ott.assert_almost_equal(ks2.computePDF(right[nDist]), expectedRight[nDist][1], 1e-5, 0.0)

    algo3 = ot.KernelSmoothing(ot.Normal(), False)
    algo3.setBoundingOption(ot.KernelSmoothing.LOWER)
    algo3.setLowerBound(baseDist.getRange().getLowerBound()[0])
    algo3.setAutomaticLowerBound(False)
    ks3 = algo3.build(sample)
    ott.assert_almost_equal(ks3.computePDF(left[nDist]), expectedLeft[nDist][2], 1e-5, 0.0)
    ott.assert_almost_equal(ks3.computePDF(right[nDist]), expectedRight[nDist][2], 1e-5, 0.0)

    algo4 = ot.KernelSmoothing(ot.Normal(), False)
    algo4.setBoundingOption(ot.KernelSmoothing.UPPER)
    algo4.setAutomaticUpperBound(True)
    ks4 = algo4.build(sample)
    ott.assert_almost_equal(ks4.computePDF(left[nDist]), expectedLeft[nDist][3], 1e-5, 0.0)
    ott.assert_almost_equal(ks4.computePDF(right[nDist]), expectedRight[nDist][3], 1e-5, 0.0)

    algo5 = ot.KernelSmoothing(ot.Normal(), False)
    algo5.setBoundingOption(ot.KernelSmoothing.UPPER)
    algo5.setUpperBound(baseDist.getRange().getUpperBound()[0])
    algo5.setAutomaticLowerBound(False)
    ks5 = algo5.build(sample)
    ott.assert_almost_equal(ks5.computePDF(left[nDist]), expectedLeft[nDist][4], 1e-5, 0.0)
    ott.assert_almost_equal(ks5.computePDF(right[nDist]), expectedRight[nDist][4], 1e-5, 0.0)

    algo6 = ot.KernelSmoothing(ot.Normal(), False)
    algo6.setBoundingOption(ot.KernelSmoothing.BOTH)
    ks6 = algo6.build(sample)
    ott.assert_almost_equal(ks6.computePDF(left[nDist]), expectedLeft[nDist][5], 1e-5, 0.0)
    ott.assert_almost_equal(ks6.computePDF(right[nDist]), expectedRight[nDist][5], 1e-5, 0.0)

    algo7 = ot.KernelSmoothing(ot.Normal(), False)
    algo7.setBoundingOption(ot.KernelSmoothing.BOTH)
    algo7.setLowerBound(baseDist.getRange().getLowerBound()[0])
    ks7 = algo7.build(sample)
    ott.assert_almost_equal(ks7.computePDF(left[nDist]), expectedLeft[nDist][6], 1e-5, 0.0)
    ott.assert_almost_equal(ks7.computePDF(right[nDist]), expectedRight[nDist][6], 1e-5, 0.0)

    algo8 = ot.KernelSmoothing(ot.Normal(), False)
    algo8.setBoundingOption(ot.KernelSmoothing.BOTH)
    algo8.setUpperBound(baseDist.getRange().getUpperBound()[0])
    ks8 = algo8.build(sample)
    ott.assert_almost_equal(ks8.computePDF(left[nDist]), expectedLeft[nDist][7], 1e-5, 0.0)
    ott.assert_almost_equal(ks8.computePDF(right[nDist]), expectedRight[nDist][7], 1e-5, 0.0)

    algo9 = ot.KernelSmoothing(ot.Normal(), False)
    algo9.setBoundingOption(ot.KernelSmoothing.BOTH)
    algo9.setLowerBound(baseDist.getRange().getLowerBound()[0])
    algo9.setUpperBound(baseDist.getRange().getUpperBound()[0])
    ks9 = algo9.build(sample)
    ott.assert_almost_equal(ks9.computePDF(left[nDist]), expectedLeft[nDist][8], 1e-5, 0.0)
    ott.assert_almost_equal(ks9.computePDF(right[nDist]), expectedRight[nDist][8], 1e-5, 0.0)

# full degenerate case
sample = ot.JointDistribution([ot.Dirac(-7.0), ot.Dirac(0.0), ot.Dirac(8.0)]).getSample(
    50
)
smoothed = ot.KernelSmoothing().build(sample)
ott.assert_almost_equal(
    smoothed.getSample(3), [[-7.0, 0.0, 8.0], [-7.0, 0.0, 8.0], [-7.0, 0.0, 8.0]], 0.0, 0.0
)

# n-d degenerate case
sample = ot.JointDistribution(
    [ot.Dirac(-7.0), ot.Arcsine(2.0, 3.0), ot.Dirac(8.0)]
).getSample(50)
sample.setDescription(["d7", "a23", "d8"])
smoothed = ot.KernelSmoothing().build(sample)
degenerateSample = smoothed.getSample(3)
ott.assert_almost_equal(
    degenerateSample, [[-7.0, 2.82969, 8.0], [-7.0, 2.35796, 8.0], [-7.0, 2.08774, 8.0]], 1e-5, 0.0
)
assert list(degenerateSample.getDescription()) == ["d7", "a23", "d8"]

# Test with reduced Cutoff - generates non positive phiGammaH
distribution = ot.Normal()
kernel = ot.Normal()
factory = ot.KernelSmoothing(kernel)
ot.ResourceMap.SetAsScalar("KernelSmoothing-CutOffPlugin", 3.0)
ot.RandomGenerator.SetSeed(8457)
sample = distribution.getSample(30)
h = factory.computePluginBandwidth(sample)[0]
ott.assert_almost_equal(h, 0.769339, 1e-5, 0.0)

# test of logTransform
for i, distribution in enumerate(
    [
        ot.LogNormal(0.0, 2.5),
        ot.Beta(20000.5, 2.5, 0.0, 1.0),
        ot.Exponential(),
        ot.WeibullMax(1.0, 0.9, 0.0),
        ot.Mixture(
            [ot.LogNormal(-1.0, 1.0, -1.0), ot.LogNormal(1.0, 1.0, 1.0)], [0.2, 0.8]
        ),
    ]
):
    sample = distribution.getSample(10000)
    kernel = ot.KernelSmoothing()
    kernel.setUseLogTransform(True)
    fitted = kernel.build(sample)
    quantile = distribution.computeQuantile(0.25)
    ott.assert_almost_equal(
        distribution.computePDF(quantile), fitted.computePDF(quantile), 0.05
    )
    quantile = distribution.computeQuantile(0.5)
    ott.assert_almost_equal(
        distribution.computePDF(quantile), fitted.computePDF(quantile), 0.05
    )
    quantile = distribution.computeQuantile(0.75)
    ott.assert_almost_equal(
        distribution.computePDF(quantile), fitted.computePDF(quantile), 0.05
    )

# segfault non-regression
data = ot.Sample.ImportFromCSVFile("t_KernelSmoothing_segfault.csv", ",")
factory = ot.KernelSmoothing()
bandwidth = [4.0 * 4.3461e09]
distribution = factory.build(data, bandwidth)

# weighted kernel smoothing, see issue #1554
ks = ot.KernelSmoothing()
twoPoints = ot.Sample([[0.0], [10.0]])
fitted = ks.buildWeighted(twoPoints, [1.0, 9.0])
ott.assert_almost_equal(fitted.getMean()[0], 9.0, 1e-10, 0.0)
# swapping the weights moves the mean accordingly
fittedSwapped = ks.buildWeighted(twoPoints, [9.0, 1.0])
ott.assert_almost_equal(fittedSwapped.getMean()[0], 1.0, 1e-10, 0.0)
# the mean is the weighted mean and the total mass is preserved
sample = ot.Normal().getSample(100)
weights = ot.Point([float(i) for i in range(100)])
totalWeight = sum(range(100))
weightedMean = (
    sum(float(i) * sample[i, 0] for i in range(100)) / totalWeight
)
fittedFull = ks.buildWeighted(sample, weights)
ott.assert_almost_equal(fittedFull.getMean()[0], weightedMean, 1e-9, 0.0)
ott.assert_almost_equal(fittedFull.computeCDF(1e8), 1.0, 1e-6)
ott.assert_almost_equal(fittedFull.computeCDF(-1e8), 0.0, 1e-6)
# zero weights are allowed
fittedZero = ks.buildWeighted(twoPoints, [0.0, 1.0])
ott.assert_almost_equal(fittedZero.getMean()[0], 10.0, 1e-10, 0.0)
# invalid calls raise
for bad in ([1.0], [-1.0, 2.0], [0.0, 0.0]):
    with ott.assert_raises(TypeError):
        ks.buildWeighted(twoPoints, bad)

# buildWeighted applies the configured log transform and binning, see #1554
sample = ot.LogNormal(2.0, 0.5).getSample(500)
equalWeights = [1.0 / 500.0] * 500
ks.setUseLogTransform(True)
fittedLog = ks.buildWeighted(sample, equalWeights)
ks.setUseLogTransform(False)
fittedPlain = ks.build(sample)
# the log-transform fit matches the plain one on the sample quantiles
for prob in (0.25, 0.5, 0.75):
    quantile = sample.computeQuantile(prob)[0]
    ott.assert_almost_equal(
        fittedLog.computeCDF([quantile]), fittedPlain.computeCDF([quantile]), 1e-3, 2e-2
    )

# invalid calls raise
with ott.assert_raises(TypeError):
    ot.KernelSmoothing(ot.Normal(2))
with ott.assert_raises(TypeError):
    ot.KernelSmoothing(ot.Normal(), True, 1)
ot.ResourceMap.SetAsUnsignedInteger("KernelSmoothing-BinNumber", 1)
with ott.assert_raises(TypeError):
    ot.KernelSmoothing()
ot.ResourceMap.SetAsUnsignedInteger("KernelSmoothing-BinNumber", 1024)

# bandwidth helpers on invalid samples
sample2d = ot.Normal(2).getSample(20)
with ott.assert_raises(TypeError):
    ot.KernelSmoothing().computePluginBandwidth(sample2d)
with ott.assert_raises(TypeError):
    ot.KernelSmoothing().computeMixedBandwidth(sample2d)
with ott.assert_raises(RuntimeError):
    ot.KernelSmoothing().computePluginBandwidth(ot.Sample(20, [1.0]))

# bandwidth/sample dimension mismatch
with ott.assert_raises(ValueError):
    ot.KernelSmoothing().build(sample2d, [0.1])

# binning is not available in dimension > 2
sample3d = ot.Normal(3).getSample(20)
with ott.assert_raises(RuntimeError):
    ot.KernelSmoothing().buildWeightedAsMixture(
        sample3d, [1.0] * 20, [0.1] * 3
    )

# boundary correction on multidimensional samples
ks2d = ot.KernelSmoothing(ot.Normal(), False)
ks2d.setBoundingOption(ot.KernelSmoothing.LOWER)
with ott.assert_raises(RuntimeError):
    ks2d.buildWeighted(sample2d, [1.0] * 20)

# user bound inconsistent with the sample
ksb = ot.KernelSmoothing(ot.Normal(), False)
ksb.setBoundingOption(ot.KernelSmoothing.LOWER)
ksb.setLowerBound(5.0)
ksb.setAutomaticLowerBound(False)
with ott.assert_raises(TypeError):
    ksb.build(ot.Normal().getSample(20))

# boundary correction on a constant sample
ksc = ot.KernelSmoothing(ot.Normal(), False)
ksc.setBoundingOption(ot.KernelSmoothing.LOWER)
with ott.assert_raises(RuntimeError):
    ksc.build(ot.Sample(20, [1.0]))

# bounding option round-trip, including the legacy constructor flag
ks = ot.KernelSmoothing()
assert ks.getBoundingOption() == ot.KernelSmoothing.NONE
ks.setBoundingOption(ot.KernelSmoothing.LOWER)
assert ks.getBoundingOption() == ot.KernelSmoothing.LOWER
ks.setBoundingOption(ot.KernelSmoothing.UPPER)
assert ks.getBoundingOption() == ot.KernelSmoothing.UPPER
legacy = ot.KernelSmoothing(ot.Normal(), False, 64, True)
assert legacy.getBoundingOption() == ot.KernelSmoothing.BOTH

# binning/kernel/log-transform accessors
ks = ot.KernelSmoothing()
assert ks.getBinning()
assert ks.getBinNumber() == 1024
ks.setBinning(False)
assert not ks.getBinning()
ks.setBinNumber(64)
assert ks.getBinNumber() == 64
with ott.assert_raises(TypeError):
    ks.setBinNumber(1)
assert ks.getKernel().getName() == "Normal"
assert not ks.getUseLogTransform()

# 2D binning path
ot.RandomGenerator.SetSeed(0)
sample = ot.Normal(2).getSample(5000)
binned2d = ot.KernelSmoothing(ot.Normal(), True, 32).build(sample)
ott.assert_almost_equal(binned2d.getMean(), [0.0, 0.0], 0.05, 0.05)

# constant and partially degenerate weighted builds
constant = ot.Sample(10, [7.0])
flat = ot.KernelSmoothing().buildWeighted(constant, [1.0] * 10, [0.5])
assert "Dirac" in repr(flat)
ott.assert_almost_equal(flat.getMean(), [7.0])
joint = ot.JointDistribution([ot.Dirac(-7.0), ot.Normal()])
jsample = joint.getSample(50)
jsmooth = ot.KernelSmoothing().buildWeighted(jsample, [1.0] * 50)
assert jsmooth.getDimension() == 2

# weighted log-transform on left-skewed data
sample = ot.WeibullMax(1.0, 0.9, 0.0).getSample(500)
ks = ot.KernelSmoothing()
ks.setUseLogTransform(True)
assert ks.getUseLogTransform()
fitted = ks.buildWeighted(sample, [1.0] * 500)
ott.assert_almost_equal(fitted.computeCDF(1e8), 1.0, 1e-6)
ott.assert_almost_equal(fitted.computeCDF(-1e8), 0.0, 1e-6)

# user upper bound inconsistent with the sample
ksu = ot.KernelSmoothing(ot.Normal(), False)
ksu.setBoundingOption(ot.KernelSmoothing.UPPER)
ksu.setUpperBound(-5.0)
ksu.setAutomaticUpperBound(False)
with ott.assert_raises(TypeError):
    ksu.build(ot.Normal().getSample(20))

# truncated build on a constant sample
ksc = ot.KernelSmoothing(ot.Normal(), False)
ksc.setBoundingOption(ot.KernelSmoothing.BOTH)
with ott.assert_raises(TypeError):
    ksc.buildAsTruncatedDistribution(ot.Sample(10, [1.0]), [0.5])

# binned truncated build
sample = ot.Uniform(-1.0, 1.0).getSample(500)
bins = ot.KernelSmoothing(ot.Normal(), True, 32)
bins.setBoundingOption(ot.KernelSmoothing.BOTH)
fitted = bins.build(sample)
ott.assert_almost_equal(fitted.computeCDF(0.0), 0.5, 0.1, 0.1)

# silverman bandwidth with zero interquartile range falls back to std
ties = ot.Sample([[0.0]] * 80 + [[1.0]] * 20)
h = ot.KernelSmoothing().computeSilvermanBandwidth(ties)
assert h[0] > 0.0
ott.assert_almost_equal(ot.KernelSmoothing().build(ties).getMean(), [0.2])

# dimension mismatch in direct builders
with ott.assert_raises(ValueError):
    ot.KernelSmoothing().buildAsKernelMixture(sample2d, [0.1])
with ott.assert_raises(ValueError):
    ot.KernelSmoothing().buildAsMixture(sample2d, [0.1])

# save/load round-trip
study = ot.Study()
study.setStorageManager(ot.XMLStorageManager("ks.xml"))
factory = ot.KernelSmoothing(ot.Normal(), True, 64)
factory.setBoundingOption(ot.KernelSmoothing.BOTH)
factory.setUseLogTransform(True)
study.add("factory", factory)
study.save()
loaded = ot.KernelSmoothing()
reloader = ot.Study()
reloader.setStorageManager(ot.XMLStorageManager("ks.xml"))
reloader.load()
reloader.fillObject("factory", loaded)
assert loaded.getBoundingOption() == ot.KernelSmoothing.BOTH
assert loaded.getUseLogTransform()
assert loaded.getBinNumber() == 64
os.remove("ks.xml")

# Non-finite bandwidths are rejected on non-degenerate marginals
ot.RandomGenerator.SetSeed(0)
bandwidthSample = ot.Normal().getSample(20)
with ott.assert_raises(TypeError):
    ot.KernelSmoothing().build(bandwidthSample, [float("inf")])
with ott.assert_raises(TypeError):
    ot.KernelSmoothing().buildWeighted(
        bandwidthSample, [1.0] * 20, [float("inf")]
    )
