"""Point 5b: cost of the discrete/continuous pairing simplification.

When both continuous and discrete atoms exist in dimension 1, each
(continuous, discrete) pair is replaced by a Mixture whose atoms are
translated copies of the continuous atom: the number of components is the
size of the discrete support. Compare with the unfused evaluation.
"""

import time

import openturns as ot

ot.TESTPREAMBLE()

atoms = [ot.Logistic(), ot.Binomial(2, 0.5), ot.Uniform()]
weights = [1.0, 2.0, 3.0]
constant = 2.0


def bench(distribution, label):
    mean = distribution.getMean()
    sigma = distribution.getStandardDeviation()
    xs = [mean[0] + (i - 2) * sigma[0] for i in range(5)]
    t0 = time.time()
    for _ in range(3):
        for x in xs:
            distribution.computePDF(x)
    t_pdf = (time.time() - t0) / 15.0
    t0 = time.time()
    for x in xs:
        distribution.computeCDF(x)
    t_cdf = (time.time() - t0) / 5.0
    coll = distribution.getDistributionCollection()
    print(
        "%-24s atoms=%2d pdf=%7.1f ms/call cdf=%7.1f ms/call"
        % (
            label,
            len(coll),
            1000.0 * t_pdf,
            1000.0 * t_cdf,
        )
    )


ot.ResourceMap.SetAsBool("LinearCombinationDistribution-SimplifyAtoms", True)
d_fused = ot.LinearCombinationDistribution(atoms, weights, constant)
bench(d_fused, "pairing ON")
coll = d_fused.getDistributionCollection()
print("  classes:", [a.getImplementation().getClassName() for a in coll])

ot.ResourceMap.SetAsBool("LinearCombinationDistribution-SimplifyAtoms", False)
d_raw = ot.LinearCombinationDistribution(atoms, weights, constant)
bench(d_raw, "pairing OFF")

# heavier discrete atom
ot.ResourceMap.SetAsUnsignedInteger(
    "LinearCombinationDistribution-MaximumSupportSize", 10000
)
ot.ResourceMap.SetAsBool("LinearCombinationDistribution-SimplifyAtoms", True)
atoms2 = [ot.Logistic(), ot.Binomial(20, 0.5), ot.Uniform()]
d_heavy = ot.LinearCombinationDistribution(atoms2, weights, constant)
bench(d_heavy, "pairing ON (n=20)")
