#! /usr/bin/env python

import pickle
from io import BytesIO

import openturns as ot
import openturns.experimental as otexp
import openturns.testing as ott

ot.TESTPREAMBLE()
ot.RandomGenerator.SetSeed(0)

factory = ot.OrthogonalProductPolynomialFactory([ot.Uniform(-1.0, 1.0)])
basis = ot.OrthogonalBasis(factory)

# sizing rule: largest m with 10*m*log(m) <= n
assert otexp.SubsampleChristoffelExperiment.DeduceSpaceDimension(0) == 1
assert otexp.SubsampleChristoffelExperiment.DeduceSpaceDimension(10) == 1
assert otexp.SubsampleChristoffelExperiment.DeduceSpaceDimension(60) == 4
assert otexp.SubsampleChristoffelExperiment.DeduceSpaceDimension(100) == 5

experiment = otexp.SubsampleChristoffelExperiment(basis, 60)
assert experiment.getSize() == 60
assert experiment.getSpaceDimension() == 4
assert experiment.getDistribution().getDimension() == 1
assert not experiment.hasUniformWeights()
assert experiment.isRandom()
assert "SubsampleChristoffelExperiment" in repr(experiment)

# weighted sample: n points in the range, positive weights of mean about 1
sample, weights = experiment.generateWithWeights()
assert sample.getSize() == 60
assert len(weights) == 60
assert all(w > 0.0 for w in weights)
assert abs(sum(weights) - 60.0) < 3.0
for point in sample:
    assert -1.0 <= point[0] <= 1.0

# default Removal: weights are exactly the density ratios, no reweighting
reference = experiment.getDistribution()
christoffel = otexp.ChristoffelDistribution(basis, 4)
ott.assert_almost_equal(weights[0], reference.computePDF(sample[0]) / christoffel.computePDF(sample[0]))

# the thinned design satisfies the frame bounds
basisFunctions = [basis.build(j) for j in range(4)]
gram = ot.SymmetricMatrix(4)
for a in range(4):
    for b in range(a + 1):
        gram[a, b] = sum(weights[k] * basisFunctions[a](sample[k])[0] * basisFunctions[b](sample[k])[0] for k in range(60)) / 60.0
for eigenvalue in gram.computeEigenValues():
    assert 0.5 <= eigenvalue <= 1.5

# setters rebuild the derived members
experiment.setSize(100)
assert experiment.getSize() == 100
assert experiment.getSpaceDimension() == 5
experiment.setBasis(basis)
assert experiment.getSpaceDimension() == 5

# the thinning is deterministic: same seed gives the same design
ot.RandomGenerator.SetSeed(0)
sampleA, weightsA = experiment.generateWithWeights()
ot.RandomGenerator.SetSeed(0)
sampleB, weightsB = experiment.generateWithWeights()
ott.assert_almost_equal(sampleA, sampleB)
ott.assert_almost_equal(weightsA, weightsB)

# the sampling factor key drives the sizing rule
ot.ResourceMap.SetAsScalar("SubsampleChristoffelExperiment-SamplingFactor", 5.0)
try:
    assert otexp.SubsampleChristoffelExperiment.DeduceSpaceDimension(60) == 6
finally:
    ot.ResourceMap.SetAsScalar("SubsampleChristoffelExperiment-SamplingFactor", 10.0)
ott.assert_almost_equal(ot.ResourceMap.GetAsScalar("SubsampleChristoffelExperiment-SamplingFactor"), 10.0)
ott.assert_almost_equal(ot.ResourceMap.GetAsScalar("SubsampleChristoffelExperiment-PoolOversamplingFactor"), 2.0)
ott.assert_almost_equal(ot.ResourceMap.GetAsScalar("SubsampleChristoffelExperiment-FrameTolerance"), 0.5)

# Barrier method: non-trivial reweighting, smoke only (no frame bounds)
ot.ResourceMap.SetAsString("SubsampleChristoffelExperiment-ThinningMethod", "Barrier")
try:
    barrierSample, barrierWeights = experiment.generateWithWeights()
    assert barrierSample.getSize() == 100
    assert len(barrierWeights) == 100
    assert all(w > 0.0 for w in barrierWeights)
finally:
    ot.ResourceMap.SetAsString("SubsampleChristoffelExperiment-ThinningMethod", "Removal")

# small space dimension (m=2): barrier schedule stays inside the spectrum
smallExperiment = otexp.SubsampleChristoffelExperiment(basis, 30)
assert smallExperiment.getSpaceDimension() == 2
smallSample, smallWeights = smallExperiment.generateWithWeights()
assert smallSample.getSize() == 30
assert len(smallWeights) == 30
assert all(w > 0.0 for w in smallWeights)

# pool factor 1 fast path: no thinning, direct Christoffel draws
ot.ResourceMap.SetAsScalar("SubsampleChristoffelExperiment-PoolOversamplingFactor", 1.0)
try:
    directSample, directWeights = experiment.generateWithWeights()
    assert directSample.getSize() == 100
    assert len(directWeights) == 100
finally:
    ot.ResourceMap.SetAsScalar("SubsampleChristoffelExperiment-PoolOversamplingFactor", 2.0)

# clone, comparison and persistence
assert experiment == otexp.SubsampleChristoffelExperiment(basis, 100)
assert experiment != otexp.SubsampleChristoffelExperiment(basis, 60)
buf = BytesIO()
pickle.dump(experiment, buf)
buf.seek(0)
assert pickle.load(buf) == experiment

# errors: the reference distribution is embedded in the basis
with ott.assert_raises(TypeError):
    experiment.setDistribution(ot.Normal())
with ott.assert_raises(TypeError):
    otexp.SubsampleChristoffelExperiment(basis, 0)
with ott.assert_raises(TypeError):
    experiment.setSize(0)
