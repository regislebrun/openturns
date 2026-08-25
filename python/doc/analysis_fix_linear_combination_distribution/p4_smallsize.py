"""Point 4: the SmallSize threshold only applies to computeProbability.

For 1D mixtures whose collection size reaches
LinearCombinationDistribution-SmallSize, computeProbability switches to the
generic algorithm while computePDF/computeCDF keep using Poisson's
summation formula regardless of the collection size.
"""

import time

import numpy as np
import openturns as ot

ot.TESTPREAMBLE()
ot.ResourceMap.SetAsBool("LinearCombinationDistribution-SimplifyAtoms", False)
ot.ResourceMap.SetAsUnsignedInteger("LinearCombinationDistribution-SmallSize", 100)

size = 1000000
print("sum of N uniform U(0,1) atoms, MC reference on %d samples" % size)
for n in [99, 100, 101]:
    d = ot.LinearCombinationDistribution([ot.Uniform(0.0, 1.0)] * n)
    mean = d.getMean()[0]
    sigma = d.getStandardDeviation()[0]
    x = mean + sigma
    sample = np.asarray(d.getSample(size))[:, 0]
    mc_cdf = float((sample <= x).mean())
    t0 = time.time()
    cdf = d.computeCDF(x)
    t_cdf = time.time() - t0
    t0 = time.time()
    prob = d.computeProbability(ot.Interval(mean - sigma, mean + sigma))
    t_prob = time.time() - t0
    t0 = time.time()
    pdf = d.computePDF(x)
    t_pdf = time.time() - t0
    # consistency: P([x-h, x+h]) vs pdf(x)*2h
    h = 0.01
    prob_local = d.computeProbability(ot.Interval(x - h, x + h))
    print(
        "N=%3d | cdf(%.3f)=%.6f (MC %.6f) %+.2e | P(m+-s)=%.6f %.3fs | "
        "pdf=%.5f %.3fs | pdf*2h vs local P: %.2e"
        % (
            n,
            x,
            cdf,
            mc_cdf,
            cdf - mc_cdf,
            prob,
            t_prob,
            pdf,
            t_pdf,
            abs(prob_local - pdf * 2 * h),
        )
    )
