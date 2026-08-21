#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott

ot.TESTPREAMBLE()
ot.RandomGenerator.SetSeed(0)

# First, compute the volume of the unit ball in R^n
a = -1.0
b = 1.0
formula = "1.0"
lower = list()
upper = list()
algo = ot.IteratedQuadrature(
    ot.GaussKronrod(20, 1.0e-6, ot.GaussKronrodRule(ot.GaussKronrodRule.G3K7))
)
for n in range(3):
    inVars = ot.Description.BuildDefault(n + 1, "x")
    inVarsBounds = inVars[0:n]
    if n > 0:
        formula += "-" + inVars[n - 1] + "^2"
        lower.append(ot.SymbolicFunction(inVarsBounds, ["-sqrt(" + formula + ")"]))
        upper.append(ot.SymbolicFunction(inVarsBounds, ["sqrt(" + formula + ")"]))
    integrand = ot.SymbolicFunction(inVars, ["1.0"])
    value = algo.integrate(integrand, a, b, lower, upper)[0]
    print(
        "dim=", n + 1, ", volume= %.12g" % value, ", calls=", integrand.getCallsNumber()
    )
# Second, integrate a multi-valued function
bounds = ot.Interval([-1.0] * 3, [1.0] * 3)
vars = ["x0", "x1", "x2"]
formulas = ["x0^2 + 2*x1^2 + 3*x2^2", "x2^2 + 2*x1^2 + 3*x0^2"]
integrand = ot.SymbolicFunction(vars, formulas)
value = algo.integrate(integrand, bounds)
print("value=", value, ", calls=", integrand.getCallsNumber())
print("Algo is based on", algo.getAlgorithm())
algo.setAlgorithm(ot.GaussLegendre([10]))
print("Algo is now based on", algo.getAlgorithm())

# Third, a 4D integration over nested bounds: the recursion depth is high
# enough to exercise the batch evaluation and the innermost integrand reuse
lo = [
    ot.SymbolicFunction(["x0"], ["0.0"]),
    ot.SymbolicFunction(["x0", "x1"], ["0.0"]),
    ot.SymbolicFunction(["x0", "x1", "x2"], ["0.0"]),
]
up = [
    ot.SymbolicFunction(["x0"], ["1.0"]),
    ot.SymbolicFunction(["x0", "x1"], ["1.0"]),
    ot.SymbolicFunction(["x0", "x1", "x2"], ["1.0"]),
]
integrand = ot.SymbolicFunction(["x0", "x1", "x2", "x3"], ["1.0"])
value = algo.integrate(integrand, 0.0, 1.0, lo, up)
ott.assert_almost_equal(value, [1.0])

# Fourth, non-rectangular bounds with a parallel evaluable integrand: the
# SymbolicEvaluation is parallel so the local integrals are computed on several
# threads, each with its own integrand and GaussKronrod workspace
algo = ot.IteratedQuadrature(
    ot.GaussKronrod(20, 1.0e-6, ot.GaussKronrodRule(ot.GaussKronrodRule.G3K7))
)
lo = [
    ot.SymbolicFunction(["x0"], ["-sqrt(1-x0^2)"]),
    ot.SymbolicFunction(["x0", "x1"], ["-sqrt(1-x0^2-x1^2)"]),
]
up = [
    ot.SymbolicFunction(["x0"], ["sqrt(1-x0^2)"]),
    ot.SymbolicFunction(["x0", "x1"], ["sqrt(1-x0^2-x1^2)"]),
]
integrand = ot.SymbolicFunction(["x0", "x1", "x2"], ["x0*x1*x2"])
assert integrand.getImplementation().isParallel()
value = algo.integrate(integrand, -1.0, 1.0, lo, up)
ott.assert_almost_equal(value, [0.0])
# Same with several outputs
integrand = ot.SymbolicFunction(["x0", "x1", "x2"], ["x0*x1", "x1*x2"])
value = algo.integrate(integrand, -1.0, 1.0, lo, up)
ott.assert_almost_equal(value, [0.0, 0.0])

# Fifth, a non GaussKronrod inner algorithm, exercising the generic path of
# the parallel evaluation
algo.setAlgorithm(ot.GaussLegendre([10]))
integrand = ot.SymbolicFunction(["x0", "x1"], ["x0*x1"])
value = algo.integrate(integrand, ot.Interval([0.0] * 2, [1.0] * 2))
ott.assert_almost_equal(value, [0.25])

# Invalid inputs are rejected
with ott.assert_raises(TypeError):
    algo.integrate(
        ot.SymbolicFunction(["x0", "x1"], ["x0*x1"]),
        0.0,
        1.0,
        [ot.SymbolicFunction(["x0"], ["0.0"])],
        [],
    )
with ott.assert_raises(TypeError):
    algo.integrate(
        ot.SymbolicFunction(["x0", "x1", "x2"], ["x0*x1*x2"]),
        0.0,
        1.0,
        [ot.SymbolicFunction(["x0"], ["0.0"])],
        [ot.SymbolicFunction(["x0"], ["1.0"])],
    )


# A non-finite value produced by the integrand must be detected
def nonFiniteFunction(X):
    return [float("nan") if X[1] > 0.5 else 0.0]


with ott.assert_raises(RuntimeError):
    algo.integrate(
        ot.PythonFunction(2, 1, nonFiniteFunction),
        0.0,
        1.0,
        [ot.SymbolicFunction(["x0"], ["0.0"])],
        [ot.SymbolicFunction(["x0"], ["1.0"])],
    )


def nonFiniteFunctionInf(X):
    return [float("inf") if X[1] > 0.5 else 0.0]


with ott.assert_raises(RuntimeError):
    algo.integrate(
        ot.PythonFunction(2, 1, nonFiniteFunctionInf),
        0.0,
        1.0,
        [ot.SymbolicFunction(["x0"], ["0.0"])],
        [ot.SymbolicFunction(["x0"], ["1.0"])],
    )

# Sixth, repr/str do not crash
algo2 = ot.IteratedQuadrature()
r = repr(algo2)
s = str(algo2)
assert "IteratedQuadrature" in s

# Seventh, 1D integration collapses to the inner algorithm directly
algo3 = ot.IteratedQuadrature(
    ot.GaussKronrod(20, 1.0e-6, ot.GaussKronrodRule(ot.GaussKronrodRule.G3K7))
)
integrand1d = ot.SymbolicFunction(["x0"], ["x0^2"])
value = algo3.integrate(integrand1d, 0.0, 1.0, [], [])
ott.assert_almost_equal(value, [1.0 / 3.0])
# Same via the Interval overload
value2 = algo3.integrate(integrand1d, ot.Interval([0.0], [1.0]))
ott.assert_almost_equal(value2, [1.0 / 3.0])

# Eighth, non-constant rectangular bounds via 3D interval
algo4 = ot.IteratedQuadrature(
    ot.GaussKronrod(20, 1.0e-6, ot.GaussKronrodRule(ot.GaussKronrodRule.G3K7))
)
integrand3d = ot.SymbolicFunction(["x0", "x1", "x2"], ["1.0"])
value = algo4.integrate(integrand3d, ot.Interval([0.0] * 3, [2.0] * 3))
ott.assert_almost_equal(value, [8.0])
