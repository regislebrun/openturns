#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import math

ot.TESTPREAMBLE()

# First, a smooth function
f = ot.SymbolicFunction("x", "sin(x)")
a = -2.5
b = 4.5
# Default parameters
algo = ot.GaussKronrod()
rules = ot.GaussKronrod.GetRules()
for i in range(len(rules)):
    algo.setRule(rules[i])
    print("Algo=", algo)
    # High-level interface
    error = ot.Point()
    value = algo.integrate(f, ot.Interval(a, b), error)[0]
    ref = math.cos(a) - math.cos(b)
    print(
        "value=%.6f" % value,
        ", ref=%.6f" % ref,
        ", true error below bound? ",
        abs(ref - value) < algo.getMaximumError(),
        ", estimated error below bound? ",
        error[0] < algo.getMaximumError(),
    )

# Second, a piecewise smooth function
f = ot.SymbolicFunction("x", "abs(sin(x))")
a = -2.5
b = 4.5
algo = ot.GaussKronrod()
rules = ot.GaussKronrod.GetRules()
rules[0] = ot.GaussKronrod.GetRuleFromName("G3K7")

for i in range(len(rules)):
    algo.setRule(rules[i])
    print("Algo=", algo)
    error = ot.Point()
    value = algo.integrate(f, ot.Interval(a, b), error)[0]
    ref = 4.0 + math.cos(b) - math.cos(a)
    print(
        "value=%.6f" % value,
        ", ref=%.6f" % ref,
        ", true error below bound? ",
        abs(ref - value) < algo.getMaximumError(),
        ", estimated error below bound? ",
        error[0] < algo.getMaximumError(),
    )

# Third, the low-level interface: containers are filled on output
algo = ot.GaussKronrod(100, 1e-12, ot.GaussKronrodRule(ot.GaussKronrodRule.G3K7))
f = ot.SymbolicFunction("x", "sin(x)")
ref = math.cos(-2.5) - math.cos(4.5)
error = ot.Point()
ai = ot.Point()
bi = ot.Point()
fi = ot.Sample()
ei = ot.Point()
value = algo.integrate(f, -2.5, 4.5, error, ai, bi, fi, ei)
ott.assert_almost_equal(value[0], ref, 1e-8, 1e-12)
ott.assert_almost_equal(error[0], 0.0, 0.0, algo.getMaximumError())
# All containers must be populated
assert ai.getDimension() > 0
assert bi.getDimension() > 0
assert fi.getSize() > 0
assert ei.getDimension() > 0

# Reusing the same containers on a second interval must not corrupt them
value = algo.integrate(f, 0.0, 1.0, error, ai, bi, fi, ei)
ref = 1.0 - math.cos(1.0)
ott.assert_almost_equal(value[0], ref, 1e-8, 1e-12)
ott.assert_almost_equal(error[0], 0.0, 0.0, algo.getMaximumError())

# Fourth, multi-output integrand through the low-level interface
g = ot.SymbolicFunction(["x"], ["sin(x)", "cos(x)"])
value2 = algo.integrate(g, -2.5, 4.5, error, ai, bi, fi, ei)
ref2 = [math.cos(-2.5) - math.cos(4.5), math.sin(4.5) - math.sin(-2.5)]
ott.assert_almost_equal(value2, ref2, 1e-8, 1e-12)
ott.assert_almost_equal(error[0], 0.0, 0.0, algo.getMaximumError())

# Fifth, Python function through the high-level interface
f_py = ot.PythonFunction(1, 1, lambda x: [math.sin(x[0])])
value3 = algo.integrate(f_py, ot.Interval(-2.5, 4.5))
ref = math.cos(-2.5) - math.cos(4.5)
ott.assert_almost_equal(value3[0], ref, 1e-8, 1e-12)

# Sixth, accessor tests
algo2 = ot.GaussKronrod(200, 1e-10, ot.GaussKronrodRule(ot.GaussKronrodRule.G7K15))
assert algo2.getMaximumSubIntervals() == 200
ott.assert_almost_equal(algo2.getMaximumError(), 1e-10)
assert "G7K15" in str(algo2.getRule())
algo2.setMaximumSubIntervals(50)
assert algo2.getMaximumSubIntervals() == 50
algo2.setMaximumError(1e-8)
ott.assert_almost_equal(algo2.getMaximumError(), 1e-8)
algo2.setRule(ot.GaussKronrodRule(ot.GaussKronrodRule.G15K31))
assert "G15K31" in str(algo2.getRule())

# Seventh, GetRules and GetRuleFromName
rules = ot.GaussKronrod.GetRules()
assert len(rules) == 6
for r in rules:
    assert r.getOrder() > 0
name_rule = ot.GaussKronrod.GetRuleFromName("G3K7")
assert "G3K7" in str(name_rule)

# Eighth, invalid inputs are rejected
with ott.assert_raises(TypeError):
    ot.GaussKronrod(1, 1e-6, ot.GaussKronrodRule(ot.GaussKronrodRule.G3K7))
with ott.assert_raises(TypeError):
    ot.GaussKronrod.GetRuleFromName("UnknownRule")
