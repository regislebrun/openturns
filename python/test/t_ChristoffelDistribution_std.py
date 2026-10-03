#! /usr/bin/env python

import math
import pickle
from io import BytesIO

import openturns as ot
import openturns.experimental as otexp
import openturns.testing as ott

ot.TESTPREAMBLE()
ot.RandomGenerator.SetSeed(0)

# default constructor is a valid 1D state (needed by doc plots)
defaultDistribution = otexp.ChristoffelDistribution()
assert defaultDistribution.getDimension() == 1
assert defaultDistribution.getSize() == 1
ott.assert_almost_equal(defaultDistribution.computePDF([0.0]), 0.5)
assert isinstance(defaultDistribution.drawPDF(), ot.Graph)
assert isinstance(defaultDistribution.drawCDF(), ot.Graph)

# orthonormal Legendre basis on Uniform(-1, 1)
factory = ot.OrthogonalProductPolynomialFactory([ot.Uniform(-1.0, 1.0)])
basis = ot.OrthogonalBasis(factory)
distribution = otexp.ChristoffelDistribution(basis, 3)

# accessors
assert distribution.getDimension() == 1
assert distribution.getSize() == 3
assert distribution.getMeasure().getDimension() == 1
assert distribution.getBasis().getSize() == 3
assert distribution.isContinuous()
assert not distribution.isDiscrete()
assert "ChristoffelDistribution" in repr(distribution)

# Christoffel function: k_3(0.5) = 1 + 3*0.5**2 + 5*P_2(0.5)**2
# with P_2(0.5) = -0.125 (orthonormal Legendre: L_2 = sqrt(5)*P_2),
# hence k_3(0.5) = 1.828125
ott.assert_almost_equal(distribution.computeChristoffel([0.5]), 1.828125)
ott.assert_almost_equal(distribution.computeChristoffel(ot.Sample([[0.5]]))[0, 0], 1.828125)

# PDF: 0.5 * 1.828125 / 3
ott.assert_almost_equal(distribution.computePDF([0.5]), 0.3046875)
ott.assert_almost_equal(distribution.computeLogPDF([0.5]), math.log(0.3046875))
assert distribution.computePDF([2.0]) == 0.0

# stability factor: 1 + 3 + 5 = 9.0 at the endpoints, Monte-Carlo estimate from below
kn = distribution.computeKn()
assert 8.0 < kn <= 9.0

# the PDF integrates to 1 over the range
probability = ot.GaussLegendre([41]).integrate(distribution.getPDF(), ot.Interval([-1.0], [1.0]))[0]
ott.assert_almost_equal(probability, 1.0)

# sampling stays inside the range
sample = distribution.getSample(10)
assert sample.getSize() == 10
assert sample.getDimension() == 1
for point in sample:
    assert -1.0 <= point[0] <= 1.0

# sequential draws target a symmetric law of mean 0
bigSample = distribution.getSample(200)
assert abs(bigSample.computeMean()[0]) < 0.2

# clone, comparison and persistence
assert distribution == otexp.ChristoffelDistribution(basis, 3)
assert distribution != otexp.ChristoffelDistribution(basis, 2)
buf = BytesIO()
pickle.dump(distribution, buf)
buf.seek(0)
ott.assert_almost_equal(pickle.load(buf).computePDF([0.5]), 0.3046875)

# errors: null size, discrete reference, wrong point dimension
with ott.assert_raises(TypeError):
    otexp.ChristoffelDistribution(basis, 0)
poissonBasis = ot.OrthogonalBasis(ot.OrthogonalProductPolynomialFactory([ot.Poisson(2.0)]))
with ott.assert_raises(TypeError):
    otexp.ChristoffelDistribution(poissonBasis, 2)
with ott.assert_raises(TypeError):
    distribution.computePDF([0.5, 0.5])
with ott.assert_raises(RuntimeError):
    distribution.getParametersCollection()

# documented ResourceMap defaults
assert ot.ResourceMap.GetAsString("ChristoffelDistribution-OptimizationAlgorithm") == "Cobyla"
assert ot.ResourceMap.GetAsUnsignedInteger("ChristoffelDistribution-RatioUniformCandidateNumber") == 10000
assert ot.ResourceMap.GetAsUnsignedInteger("ChristoffelDistribution-RatioUniformMaxDimension") == 5
assert ot.ResourceMap.GetAsUnsignedInteger("ChristoffelDistribution-KnSamplingSize") == 100000
ott.assert_almost_equal(ot.ResourceMap.GetAsScalar("ChristoffelDistribution-KnSafetyFactor"), 1.0)
assert ot.ResourceMap.GetAsUnsignedInteger("ChristoffelDistribution-SliceGridSize") == 1000

# 2D tensor case: k_2(x) = 1 + 3*x0**2 with the linear enumerate
basis2 = ot.OrthogonalBasis(ot.OrthogonalProductPolynomialFactory([ot.Uniform(-1.0, 1.0), ot.Uniform(-1.0, 1.0)]))
distribution2 = otexp.ChristoffelDistribution(basis2, 2)
ott.assert_almost_equal(distribution2.computePDF([0.5, -0.3]), 0.21875)
# X1 | X0 is uniform here since k_2 ignores x1
ott.assert_almost_equal(distribution2.computeConditionalPDF(0.3, [0.2]), 0.5)
# sequential versions fan out to the scalar ones
point2 = [0.5, -0.3]
seqPDF = distribution2.computeSequentialConditionalPDF(point2)
ott.assert_almost_equal(seqPDF[0], 0.4375)
ott.assert_almost_equal(seqPDF[1], distribution2.computeConditionalPDF(point2[1], [point2[0]]))
seqCDF = distribution2.computeSequentialConditionalCDF(point2)
ott.assert_almost_equal(seqCDF[1], distribution2.computeConditionalCDF(point2[1], [point2[0]]))
# quantile inverts the conditional CDF
level = distribution2.computeConditionalCDF(0.3, [0.2])
ott.assert_almost_equal(distribution2.computeConditionalQuantile(level, [0.2]), 0.3)
# the Christoffel law itself is never independent
assert not distribution2.hasIndependentCopula()
sample2 = distribution2.getSample(5)
assert sample2.getSize() == 5
for point in sample2:
    assert -1.0 <= point[0] <= 1.0 and -1.0 <= point[1] <= 1.0

# correlated reference: ratio-of-uniforms branch (d=2 <= 5)
corr = ot.CorrelationMatrix(2)
corr[0, 1] = 0.5
measureC = ot.Normal([0.0, 0.0], [1.0, 1.0], corr)
hermite = ot.OrthogonalBasis(ot.OrthogonalProductPolynomialFactory([ot.Normal(), ot.Normal()]))
whiten = ot.SymbolicFunction(['x0', 'x1'], ['x0', '(x1 - 0.5 * x0) / 0.8660254037844386'])
corrBasis = ot.OrthogonalBasis(otexp.FiniteOrthogonalFunctionFactory([ot.ComposedFunction(hermite.build(j), whiten) for j in range(3)], measureC))
assert not corrBasis.getMeasure().hasIndependentCopula()
distributionC = otexp.ChristoffelDistribution(corrBasis, 3)
assert distributionC.computePDF([0.1, -0.2]) >= 0.0
sampleC = distributionC.getSample(5)
assert sampleC.getSize() == 5
# trace identity: mean of k_3 under mu is 3 by orthonormality
knVals = distributionC.computeChristoffel(measureC.getSample(2000))
assert 2.5 < knVals.computeMean()[0] < 3.5

# forcing the crude rejection branch through the dimension key
ot.ResourceMap.SetAsUnsignedInteger("ChristoffelDistribution-RatioUniformMaxDimension", 0)
try:
    rejSample = distribution.getSample(5)
    assert rejSample.getSize() == 5
    for point in rejSample:
        assert -1.0 <= point[0] <= 1.0
finally:
    ot.ResourceMap.SetAsUnsignedInteger("ChristoffelDistribution-RatioUniformMaxDimension", 5)
