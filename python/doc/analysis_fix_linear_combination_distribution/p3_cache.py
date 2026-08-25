"""Point 3 (b): cache economics on the reference heavy mixture.

4-atom 3D mixture (Normal, Mixture of two Normals, two Uniforms): the
pointwise PDF stops on the MaximumPDFLevel bound after tens of levels,
beyond the default cache horizon (~19 levels in dimension 3).
"""

import time

import openturns as ot

ot.TESTPREAMBLE()
ot.ResourceMap.SetAsUnsignedInteger("LinearCombinationDistribution-DefaultMaxSize", 65536)

mixture = ot.Mixture([ot.Normal(2.0, 1.0), ot.Normal(-2.0, 1.0)])
collection = ot.DistributionCollection(0)
collection.add(ot.Normal(0.0, 1.0))
collection.add(mixture)
collection.add(ot.Uniform(0, 1))
collection.add(ot.Uniform(0, 1))
W = ot.Matrix([[1.0, -0.05, 1.0, -0.5], [0.5, 1.0, -0.05, 0.3], [-0.5, -0.1, 1.2, -0.8]])


def bench(maxSize, nPoints):
    ot.ResourceMap.SetAsUnsignedInteger(
        "LinearCombinationDistribution-DefaultMaxSize", maxSize
    )
    d = ot.LinearCombinationDistribution(collection, W)
    mean = d.getMean()
    sigma = d.getStandardDeviation()
    pts = [list(mean + (0.5 + 0.04 * i) * sigma) for i in range(nPoints)]
    t0 = time.time()
    values = [d.computePDF(p) for p in pts]
    return time.time() - t0, values


for maxSize in [65536, 16777216]:
    dt5, _ = bench(maxSize, 5)
    dt20, _ = bench(maxSize, 20)
    print(
        "maxSize=%9d  %2d pts: %7.2fs (%6.0f ms/call)   %2d pts: %7.2fs (%6.0f ms/call)"
        % (maxSize, 5, dt5, 1000.0 * dt5 / 5, 20, dt20, 1000.0 * dt20 / 20)
    )
