#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import math as m

v1 = [0.0, 0.0, 0.0]
v2 = [1.0, 1.0, 1.0]
v3 = [0.0, 1.0, 1.0]
v4 = [0.0, 0.0, 1.0]
# simplex
S = [v1, v2, v3, v4]
simplicies = [[0, 1, 2, 3]]
mesh = ot.Mesh(S, simplicies)
# function to integrate
f = ot.SymbolicFunction(["x1", "x2", "x3"], ["exp(x1 + x2 + x3)"])
# integral of f
algo = ot.SimplicialCubature()
value = algo.integrate(f, mesh)[0]
print(value)
ott.assert_almost_equal(value, (m.exp(1.0) - 1.0) ** 3 / 6)

# Test with interval inferface
f = ot.SymbolicFunction(["x", "y", "z"], ["sin(x) * cos(y) * exp(z)"])
valueRef = -m.sin(1.0) * (m.cos(1.0) - 1.0) * (m.e - 1.0)
value = algo.integrate(f, ot.Interval([0.0] * 3, [1.0] * 3))
ott.assert_almost_equal(value[0], valueRef)

# Multi-output integrand over the simplex
f = ot.SymbolicFunction(["x1", "x2", "x3"], ["x1", "x2"])
value = algo.integrate(f, mesh)
ott.assert_almost_equal(value, [1.0 / 24.0, 1.0 / 12.0])

# Adaptive refinement on a peaked integrand over [0, 1]^3
f = ot.SymbolicFunction(["x1", "x2", "x3"], ["exp(-10.0 * (x1 + x2 + x3))"])
valueRef = ((1.0 - m.exp(-10.0)) / 10.0) ** 3
algo.setRule(3)
algo.setMaximumRelativeError(1e-5)
value = algo.integrate(f, ot.Interval([0.0] * 3, [1.0] * 3))
ott.assert_almost_equal(value[0], valueRef, 1e-6, 1e-12)

# Rule 4 accuracy over the simplex
f = ot.SymbolicFunction(["x1", "x2", "x3"], ["exp(x1 + x2 + x3)"])
algo.setRule(4)
algo.setMaximumRelativeError(1e-5)
value = algo.integrate(f, mesh)[0]
ott.assert_almost_equal(value, (m.exp(1.0) - 1.0) ** 3 / 6, 1e-6, 1e-12)

# Rules 1 and 2 accuracy over the simplex and over a 2D interval
f = ot.SymbolicFunction(["x1", "x2", "x3"], ["exp(x1 + x2 + x3)"])
refSimplex = (m.exp(1.0) - 1.0) ** 3 / 6
f2 = ot.SymbolicFunction(["x", "y"], ["sin(x) * cos(y)"])
refInterval2 = -m.sin(1.0) * (m.cos(1.0) - 1.0)
for rule in [1, 2]:
    algo.setRule(rule)
    algo.setMaximumRelativeError(1e-5)
    value = algo.integrate(f, mesh)[0]
    ott.assert_almost_equal(value, refSimplex, 1e-6, 1e-12)
    value = algo.integrate(f2, ot.Interval([0.0] * 2, [1.0] * 2))[0]
    ott.assert_almost_equal(value, refInterval2, 1e-6, 1e-12)

# Adaptive refinement of a multi-output integrand: each marginal is integrated
# independently with its own error estimates
f = ot.SymbolicFunction(["x1", "x2", "x3"], ["exp(-10.0 * x1)", "exp(-10.0 * x2)", "exp(-10.0 * x3)"])
valueRef = ((1.0 - m.exp(-10.0)) / 10.0) * ot.Point([1.0, 1.0, 1.0])
algo.setRule(3)
algo.setMaximumRelativeError(1e-5)
value = algo.integrate(f, ot.Interval([0.0] * 3, [1.0] * 3))
ott.assert_almost_equal(value, valueRef, 1e-6, 1e-12)

# Regression test for the deterministic selection of the simplices to refine:
# a uniform mesh of a symmetric integrand produces many simplices with equal
# errors, so the partial sort must break the ties with the simplex index to
# select a deterministic set. The result must be identical on every call and
# accurate.
meshUniform = ot.IntervalMesher([8, 8, 8]).build(ot.Interval([0.0] * 3, [1.0] * 3))
f = ot.SymbolicFunction(["x1", "x2", "x3"], ["exp(-10.0 * (x1 + x2 + x3))"])
valueRef = ((1.0 - m.exp(-10.0)) / 10.0) ** 3
algo.setRule(3)
algo.setMaximumRelativeError(1e-7)
valueRefAdap = algo.integrate(f, meshUniform)
for _ in range(3):
    value = algo.integrate(f, meshUniform)
    ott.assert_almost_equal(value, valueRefAdap)
    ott.assert_almost_equal(value[0], valueRef, 1e-6, 1e-8)
# The refinement must actually take place for the tie-breaking to be exercised
assert f.getCallsNumber() > meshUniform.getSimplicesNumber() * 16

# The batched evaluation must give the same result for any block size
savedBlockSize = ot.ResourceMap.GetAsUnsignedInteger("SimplicialCubature-EvaluationBlockSize")
try:
    f = ot.SymbolicFunction(["x1", "x2", "x3"], ["exp(-10.0 * (x1 + x2 + x3))"])
    algo.setMaximumRelativeError(1e-6)
    reference = None
    for blockSize in [2048, 7, 1]:
        ot.ResourceMap.SetAsUnsignedInteger("SimplicialCubature-EvaluationBlockSize", blockSize)
        value = algo.integrate(f, meshUniform)
        ott.assert_almost_equal(value[0], valueRef, 1e-6, 1e-8)
        if reference is None:
            reference = value
        ott.assert_almost_equal(value, reference, 0.0, 1e-12)
finally:
    ot.ResourceMap.SetAsUnsignedInteger("SimplicialCubature-EvaluationBlockSize", savedBlockSize)

# Invalid inputs are rejected
ot.ResourceMap.SetAsUnsignedInteger("SimplicialCubature-EvaluationBlockSize", 2048)
with ott.assert_raises(TypeError):
    algo.integrate(ot.SymbolicFunction(["x"], ["x"]), mesh)
with ott.assert_raises(TypeError):
    algo.setRule(5)
with ott.assert_raises(TypeError):
    ot.ResourceMap.SetAsUnsignedInteger("SimplicialCubature-EvaluationBlockSize", 0)
    algo.integrate(f, meshUniform)
ot.ResourceMap.SetAsUnsignedInteger("SimplicialCubature-EvaluationBlockSize", 2048)

# An empty mesh of the right dimension integrates to zero
f = ot.SymbolicFunction(["x1", "x2"], ["x1", "x2"])
emptyMesh = ot.Mesh(ot.Sample(0, 2), ot.IndicesCollection([]))
value = algo.integrate(f, emptyMesh)
ott.assert_almost_equal(value, [0.0, 0.0])

# The maximum calls number bounds the number of evaluations
f = ot.SymbolicFunction(["x1", "x2", "x3"], ["exp(-20.0 * (x1 + x2 + x3))"])
algo.setRule(3)
algo.setMaximumCallsNumber(200)
algo.integrate(f, ot.Interval([0.0] * 3, [1.0] * 3))
assert f.getCallsNumber() >= 1
assert f.getCallsNumber() <= 4 * 200

# Accessor round-trip tests
algo2 = ot.SimplicialCubature()
algo2.setRule(2)
assert algo2.getRule() == 2
algo2.setMaximumAbsoluteError(1e-10)
ott.assert_almost_equal(algo2.getMaximumAbsoluteError(), 1e-10)
algo2.setMaximumRelativeError(1e-8)
ott.assert_almost_equal(algo2.getMaximumRelativeError(), 1e-8)
algo2.setMaximumCallsNumber(500)
assert algo2.getMaximumCallsNumber() == 500

# repr/str do not crash
r = repr(algo2)
s = str(algo2)
assert "SimplicialCubature" in s
assert "rule" in s

# Default constructor uses ResourceMap defaults
algo3 = ot.SimplicialCubature()
assert algo3.getRule() == ot.ResourceMap.GetAsUnsignedInteger("SimplicialCubature-DefaultRule")
ott.assert_almost_equal(algo3.getMaximumAbsoluteError(), ot.ResourceMap.GetAsScalar("SimplicialCubature-DefaultMaximumAbsoluteError"))
ott.assert_almost_equal(algo3.getMaximumRelativeError(), ot.ResourceMap.GetAsScalar("SimplicialCubature-DefaultMaximumRelativeError"))
assert algo3.getMaximumCallsNumber() == ot.ResourceMap.GetAsUnsignedInteger("SimplicialCubature-DefaultMaximumCallsNumber")

# 1D integration over an interval
f1d = ot.SymbolicFunction(["x"], ["x^2"])
value1d = algo2.integrate(f1d, ot.Interval([0.0], [1.0]))
ott.assert_almost_equal(value1d[0], 1.0 / 3.0)

# 2D integration over a triangle
f2d = ot.SymbolicFunction(["x", "y"], ["x + y"])
v1 = [0.0, 0.0]
v2 = [1.0, 0.0]
v3 = [0.0, 1.0]
mesh2d = ot.Mesh([v1, v2, v3], [[0, 1, 2]])
value2d = algo2.integrate(f2d, mesh2d)
ott.assert_almost_equal(value2d[0], 1.0 / 3.0)
