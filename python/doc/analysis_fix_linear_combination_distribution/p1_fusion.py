"""Point 1: exact-equality merging of fused discrete supports.

The fusion of discrete atoms convolves their supports scaled by their
weights and merges equal support points with an exact '==' comparison.
With several atoms, the same mathematical point can be produced through
different accumulation orders and remain duplicated when the floating
point roundings differ.
"""

from itertools import product

import openturns as ot

ot.TESTPREAMBLE()


def fused_support_size(atoms, weights):
    d = ot.LinearCombinationDistribution(atoms, weights)
    coll = d.getDistributionCollection()
    total = 0
    for a in coll:
        if a.getImplementation().getClassName() == "FiniteDiscreteDistribution":
            total += a.getSupport().getSize()
        else:
            total += a.getSupport().getSize()
    return coll, total


def exact_cardinality(atoms, weights):
    """Cardinality of the exact support of the combination."""
    values = set()
    ranges = []
    for i, a in enumerate(atoms):
        sup = [weights[i] * a.getSupport()[k, 0] for k in range(a.getSupport().getSize())]
        ranges.append(sup)
    for combo in product(*ranges):
        values.add(round(sum(combo), 9))
    return len(values)


print("--- search for duplicated points missed by the exact comparison")
cases = [
    ([ot.Binomial(2, 0.1)] * 3, [1.0, 1.0, 2.0]),
    ([ot.Binomial(2, 0.1)] * 3, [1.0, 2.0, 3.0]),
    ([ot.Binomial(2, 0.1)] * 4, [1.0, 2.0, 0.5, 0.25]),
    ([ot.Binomial(3, 0.25)] * 3, [1.0, 0.5, 0.25]),
    ([ot.Poisson(2.0)] * 3, [1.0, 1.0, 2.0]),
]
ot.ResourceMap.SetAsUnsignedInteger(
    "LinearCombinationDistribution-MaximumSupportSize", 1000000
)
for atoms, weights in cases:
    coll, size = fused_support_size(atoms, weights)
    exact = exact_cardinality(atoms, weights)
    print(
        "weights=%-22s stored=%3d exact=%3d %s"
        % (
            str(weights),
            size,
            exact,
            "(INFLATED)" if size > exact else "",
        )
    )

print()
print("--- worst case found by scanning fractional weights")
worst = None
for wi in range(1, 12):
    weights = [1.0, 0.1 * wi, 0.25 * wi]
    atoms = [ot.Binomial(2, 0.1)] * 3
    coll, size = fused_support_size(atoms, weights)
    exact = exact_cardinality(atoms, weights)
    if size > exact:
        print("INFLATED: weights=%s stored=%d exact=%d" % (weights, size, exact))
        worst = (weights, size, exact)
print("worst:", worst)

print()
print("--- impact on the fusion itself (MaximumSupportSize sweep)")
atoms = [ot.Binomial(2, 0.1)] * 10
weights = [1.0 + 0.1 * i for i in range(10)]
coll, size = fused_support_size(atoms, weights)
exact = exact_cardinality(atoms, weights)
print("10 fractional-weight binomials: stored support=%d exact=%d" % (size, exact))
for mss in [500, 1000, 5000, 100000]:
    ot.ResourceMap.SetAsUnsignedInteger(
        "LinearCombinationDistribution-MaximumSupportSize", mss
    )
    d = ot.LinearCombinationDistribution(atoms, weights)
    coll = d.getDistributionCollection()
    classes = [a.getImplementation().getClassName() for a in coll]
    nfdd = classes.count("FiniteDiscreteDistribution")
    sizes = [
        a.getSupport().getSize()
        for a in coll
        if a.getImplementation().getClassName() == "FiniteDiscreteDistribution"
    ]
    print(
        "MaximumSupportSize=%6d -> %d fused atom(s), support sizes=%s"
        % (mss, nfdd, sizes)
    )
