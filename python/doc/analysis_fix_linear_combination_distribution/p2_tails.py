"""Point 2: range truncation versus heavy-tailed atoms.

The range of a LinearCombinationDistribution is clipped to
[position indicator +/- beta * dispersion indicator] even when the exact
range is unbounded: beyond it the CDF returns exactly 1 and the PDF exactly
0, regardless of the remaining mass.
"""

import numpy as np
import openturns as ot

ot.TESTPREAMBLE()
ot.RandomGenerator.SetSeed(0)

size = 500000
print("atoms: LogNormal(muLog=0, sigmaLog) + Normal(0, 0.1), weights [1, 0.1]")
print("sigmaLog | upperBound | MC mass above | quantile(1-1e-3) == ub ?")
for sigmaLog in [1.0, 1.5, 2.0, 2.5]:
    atoms = [ot.LogNormal(0.0, sigmaLog), ot.Normal()]
    d = ot.LinearCombinationDistribution(atoms, [1.0, 0.1])
    ub = d.getRange().getUpperBound()[0]
    sample = np.asarray(d.getSample(size))
    exceedance = float((sample[:, 0] > ub).mean())
    q = d.computeQuantile(1.0 - 1.0e-3)[0]
    print(
        "%8.1f | %11.4g | %13.3g | %.6g %s"
        % (sigmaLog, ub, exceedance, q, "(clipped!)" if abs(q - ub) < 1e-12 else "")
    )

print()
print("--- same for a mixture with a Pareto-like tail (Beta(alpha<2))")
for alpha in [2.5, 1.5]:
    atoms = [ot.Pareto(alpha, 1.0), ot.Normal()]
    d = ot.LinearCombinationDistribution(atoms, [1.0, 0.1])
    ub = d.getRange().getUpperBound()[0]
    sample = np.asarray(d.getSample(size))
    exceedance = float((sample[:, 0] > ub).mean())
    q = d.computeQuantile(1.0 - 1.0e-3)[0]
    print(
        "alpha=%3.1f | upperBound=%10.4g | MC mass above=%8.3g | q=%.6g %s"
        % (alpha, ub, exceedance, q, "(clipped!)" if abs(q - ub) < 1e-12 else "")
    )
