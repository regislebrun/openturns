"""Point 5a: opposite-weight Poisson atoms could fuse into a Skellam atom.

The fusion of Poisson atoms only handles same-sign weights (poissonMap
accumulates lambda per weight). Opposite weights are left untouched although
sum_i Poi(l_i) - sum_j Poi(m_j) is a Skellam(sum l, sum m) when all the
weights are +/-1.
"""


import openturns as ot

ot.TESTPREAMBLE()
ot.ResourceMap.SetAsBool("LinearCombinationDistribution-SimplifyAtoms", True)

print("--- current behaviour")
d = ot.LinearCombinationDistribution([ot.Poisson(2.0), ot.Poisson(3.0)], [1.0, -1.0])
coll = d.getDistributionCollection()
print(
    "atoms:",
    [(a.getImplementation().getClassName(), a.getParameter()) for a in coll],
    "constant=",
    d.getConstant(),
)

print()
print("--- Skellam equivalence check (unit weights)")
l1, l2 = 2.0, 3.0
d = ot.LinearCombinationDistribution([ot.Poisson(l1), ot.Poisson(l2)], [1.0, -1.0])
sk = ot.Skellam(l1, l2)
print("LCD   mean=%.6g var=%.6g" % (d.getMean()[0], d.getCovariance()[0, 0]))
print("Skell mean=%.6g var=%.6g" % (sk.getMean()[0], sk.getCovariance()[0, 0]))
# pointwise agreement of the probability mass on a few integers
ok = True
for k in range(-3, 8):
    pk = sk.computePDF(k)
    pk_lcd = d.computePDF(k)
    ok = ok and abs(pk - pk_lcd) < 1e-12
print("pointwise PMF agreement on [-3, 8):", ok)
print(
    "isAnalytical after Skellam substitution would be True;",
    "current atoms=%d -> generic algorithms" % len(coll),
)
