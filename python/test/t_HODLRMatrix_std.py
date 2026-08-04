#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import math as m

ot.PlatformInfo.SetNumericalPrecision(3)
ot.ResourceMap.SetAsUnsignedInteger("HODLRMatrix-MinLeafSize", 4)
ot.ResourceMap.SetAsScalar("HODLRMatrix-AssemblyEpsilon", 1.0e-6)
ot.ResourceMap.SetAsScalar("HODLRMatrix-RecompressionEpsilon", 1.0e-6)


class HODLRTestAssemblyFunction:
    def __init__(self, vertices, lengthScale=0.1):
        self.vertices = vertices
        self.lengthScale = lengthScale

    def __call__(self, i, j):
        pt1 = self.vertices[i]
        pt2 = self.vertices[j]
        diff = pt1 - pt2
        return m.exp(-diff.norm() / self.lengthScale)


def build_dense_matrix(assembly_func, size):
    K = ot.Matrix(size, size)
    for i in range(size):
        for j in range(size):
            K[i, j] = assembly_func(i, j)
    return K


def make_params(epsilon=1.0e-6, leaf=4, factorization="LU"):
    p = ot.HODLRMatrixParameters()
    p.setAssemblyEpsilon(epsilon)
    p.setRecompressionEpsilon(epsilon)
    p.setMinLeafSize(leaf)
    p.setFactorizationMethod(factorization)
    return p


def build_hodlr(vertices, af, symmetric=True, params=None):
    if params is None:
        params = make_params()
    factory = ot.HODLRMatrixFactory()
    hodlr = factory.build(vertices, 1, symmetric, params)
    hodlr.assembleReal(af, 'L')
    hodlr.factorize(params.getFactorizationMethod())
    return hodlr


def setup_20():
    vertices = ot.Sample([[i * 0.1] for i in range(20)])
    af = HODLRTestAssemblyFunction(vertices, 0.1)
    hodlr = build_hodlr(vertices, af)
    K = build_dense_matrix(af, 20)
    return vertices, af, hodlr, K


# === Test 1: getNbRows / getNbColumns ===
print("=== Test 1: getNbRows / getNbColumns ===")
vertices, af, hodlr, K = setup_20()
ott.assert_almost_equal(ot.Point([hodlr.getNbRows()]), ot.Point([20]))
ott.assert_almost_equal(ot.Point([hodlr.getNbColumns()]), ot.Point([20]))
print("  PASS")

# === Test 2: copy + re-assemble ===
print("\n=== Test 2: copy ===")
hodlr_copy = hodlr.copy()
ott.assert_almost_equal(
    ot.Point([hodlr_copy.getNbRows()]),
    ot.Point([hodlr.getNbRows()])
)
# copy() drops the tree: it reports its true condition until re-assembled
with ott.assert_raises(TypeError):
    hodlr_copy.logDeterminant()
# copy() preserves dimension but not tree -- re-assemble
hodlr_copy.assembleReal(af, 'L')
hodlr_copy.factorize('LU')
b = ot.Point(20, 1.0)
ott.assert_almost_equal(hodlr_copy.solve(b), hodlr.solve(b), 1.0e-10, 1.0e-10)
print("  PASS")

# === Test 3: __repr__ / __str__ ===
print("\n=== Test 3: __repr__ / __str__ ===")
repr_str = repr(hodlr)
str_str = str(hodlr)
assert len(repr_str) > 0, "repr() should not be empty"
assert len(str_str) > 0, "str() should not be empty"
print("  PASS")

# === Test 4: getDiagonal ===
print("\n=== Test 4: getDiagonal ===")
diag = hodlr.getDiagonal()
ott.assert_almost_equal(ot.Point([diag.getSize()]), ot.Point([20]))
# exp(-||x_i-x_j||/L) at i==j is exp(0)=1
ott.assert_almost_equal(diag, ot.Point(20, 1.0), 0.15, 0.15)
print("  PASS")

# === Test 5: solve(Point) ===
print("\n=== Test 5: solve(Point) ===")
b5 = ot.Point(20, 1.0)
x_h = hodlr.solve(b5)
x_d = K.solveLinearSystem(b5)
ott.assert_almost_equal(x_h, x_d, 1.0e-4, 1.0e-4)
print("  PASS")

# === Test 6: solve(Matrix) - multi-RHS ===
print("\n=== Test 6: solve(Matrix) ===")
B = ot.Matrix(20, 3)
for i in range(20):
    B[i, 0] = 1.0
    B[i, 1] = float(i) / 20.0
    B[i, 2] = m.sin(float(i))
X_h = hodlr.solve(B)
X_d = K.solveLinearSystem(B)
for col in range(3):
    ott.assert_almost_equal(
        ot.Point([X_h[i, col] for i in range(20)]),
        ot.Point([X_d[i, col] for i in range(20)]),
        1.0e-4, 1.0e-4
    )
print("  PASS")

# === Test 7: scale ===
print("\n=== Test 7: scale ===")
hodlr_s = hodlr.copy()
hodlr_s.assembleReal(af, 'L')
hodlr_s.factorize('LU')
scale_factor = 3.5
hodlr_s.scale(scale_factor)
# scale() rebuilds the tree with a scaled evaluator and invalidates the
# factorization until factorize() is called again
b7 = ot.Point(20, 1.0)
with ott.assert_raises(TypeError):
    hodlr_s.solve(b7)
hodlr_s.factorize('LU')
# solve(alpha*A, b) = solve(A, b)/alpha
x_s = hodlr_s.solve(b7)
x_orig = hodlr.solve(b7)
ott.assert_almost_equal(x_s, x_orig / scale_factor, 1.0e-8, 1.0e-8)
# logDet(alpha*A) = n*log(alpha) + logDet(A)
ott.assert_almost_equal(
    ot.Point([hodlr_s.logDeterminant()]),
    ot.Point([hodlr.logDeterminant() + 20 * m.log(scale_factor)]),
    1.0e-8, 1.0e-8
)
# diagonal is scaled too
ott.assert_almost_equal(
    hodlr_s.getDiagonal(), hodlr.getDiagonal() * scale_factor, 1.0e-10, 1.0e-10
)
print("  PASS")

# === Test 8: addIdentity ===
print("\n=== Test 8: addIdentity ===")
hodlr_ai = hodlr.copy()
hodlr_ai.assembleReal(af, 'L')
hodlr_ai.factorize('LU')
alpha = 2.0
hodlr_ai.addIdentity(alpha)
ott.assert_almost_equal(
    hodlr_ai.getDiagonal(), hodlr.getDiagonal() + ot.Point(20, alpha), 1.0e-10, 1.0e-10
)
print("  PASS")

# === Test 9: compressionRatio / fullrkRatio ===
print("\n=== Test 9: compressionRatio / fullrkRatio ===")
cr = hodlr.compressionRatio()
fr = hodlr.fullrkRatio()
assert cr[0] > 0, "stored entries should be positive"
assert cr[0] <= cr[1], "compressed should be <= total"
assert fr[1] > 0, "rank should be positive"
print(f"  compressionRatio: stored={cr[0]}, total={cr[1]}")
print(f"  fullrkRatio: full={fr[0]}, rank={fr[1]}")
print("  PASS")

# === Test 10: HODLRMatrixParameters getters/setters ===
print("\n=== Test 10: HODLRMatrixParameters getters ===")
params = ot.HODLRMatrixParameters()
params.setAssemblyEpsilon(1.0e-3)
params.setRecompressionEpsilon(1.0e-5)
params.setMinLeafSize(16)
ott.assert_almost_equal(
    ot.Point([params.getAssemblyEpsilon()]), ot.Point([1.0e-3])
)
ott.assert_almost_equal(
    ot.Point([params.getRecompressionEpsilon()]), ot.Point([1.0e-5])
)
ott.assert_almost_equal(ot.Point([params.getMinLeafSize()]), ot.Point([16]))
print("  PASS")

# === Test 11: setFactorizationMethod / getFactorizationMethod ===
print("\n=== Test 11: factorizationMethod ===")
params2 = ot.HODLRMatrixParameters()
params2.setFactorizationMethod("LLt")
ott.assert_almost_equal(
    ot.Point([len(params2.getFactorizationMethod())]), ot.Point([3])
)
assert params2.getFactorizationMethod() == "LLt", "factorization method should be 'LLt'"
params2.setFactorizationMethod("LU")
ott.assert_almost_equal(
    ot.Point([len(params2.getFactorizationMethod())]), ot.Point([2])
)
assert params2.getFactorizationMethod() == "LU", "factorization method should be 'LU'"
print("  PASS")

# === Test 12: LLt factorization ===
print("\n=== Test 12: LLt factorization ===")
n12 = 50
vertices12 = ot.Sample([[i * 0.02] for i in range(n12)])
af12 = HODLRTestAssemblyFunction(vertices12, 0.1)
params12 = make_params(leaf=25, factorization="LLt")
params12.setMaxRank(50)
hodlr12 = build_hodlr(vertices12, af12, params=params12)
K12 = build_dense_matrix(af12, n12)
b12 = ot.Point(n12, 1.0)
x_h12 = hodlr12.solve(b12)
x_d12 = K12.solveLinearSystem(b12)
ott.assert_almost_equal(x_h12, x_d12, 1.0e-4, 1.0e-4)
print("  PASS")

# === Test 13: solveLower ===
print("\n=== Test 13: solveLower ===")
n13 = 20
vertices13 = ot.Sample([[i * 0.05] for i in range(n13)])
af13 = HODLRTestAssemblyFunction(vertices13, 0.1)
# solveLower requires the Cholesky factor: an LU factorization must be rejected
hodlr13_lu = build_hodlr(vertices13, af13)
b13 = ot.Point(n13, 1.0)
with ott.assert_raises(TypeError):
    hodlr13_lu.solveLower(b13)
# LLt: solveLower(b) = L^{-1} b, solveLower(b, True) = L^{-T} b
params13 = make_params(leaf=25, factorization="LLt")
params13.setMaxRank(50)
hodlr13 = build_hodlr(vertices13, af13, params=params13)
K13 = ot.CovarianceMatrix(n13)
for i in range(n13):
    for j in range(n13):
        K13[i, j] = af13(i, j)
L13 = K13.computeCholesky()
x_lower = hodlr13.solveLower(b13)
x_forward = L13.solveLinearSystem(b13)
ott.assert_almost_equal(x_lower, x_forward, 1.0e-4, 1.0e-4)
x_lower_t = hodlr13.solveLower(b13, True)
x_backward = L13.transpose().solveLinearSystem(b13)
ott.assert_almost_equal(x_lower_t, x_backward, 1.0e-4, 1.0e-4)
# L^{-T} L^{-1} b must reproduce the full solve C^{-1} b
x_composed = hodlr13.solveLower(hodlr13.solveLower(b13), True)
x_full = hodlr13.solve(b13)
ott.assert_almost_equal(x_composed, x_full, 1.0e-4, 1.0e-4)
print("  PASS (Point)")

B13 = ot.Matrix(n13, 2)
for i in range(n13):
    B13[i, 0] = 1.0
    B13[i, 1] = float(i)
X_lower = hodlr13.solveLower(B13)
X_forward = L13.solveLinearSystem(B13)
ott.assert_almost_equal(X_lower, X_forward, 1.0e-4, 1.0e-4)
print("  PASS (Matrix)")

# === Test 13b: gemv with LLt factorization ===
print("\n=== Test 13b: gemv with LLt ===")
# gemv must return A * x (not L * x) when the matrix is Cholesky-factorized
x13b = ot.Point([float(i) / n12 for i in range(n12)])
y13b = ot.Point(n12, 0.0)
hodlr12.gemv('N', 1.0, x13b, 0.0, y13b)
y13b_dense = ot.Point(n12, 0.0)
for i in range(n12):
    s = 0.0
    for j in range(n12):
        s += af12(i, j) * x13b[j]
    y13b_dense[i] = s
ott.assert_almost_equal(y13b, y13b_dense, 1.0e-4, 1.0e-4)
print("  PASS")

# === Test 14: solve with trans ===
print("\n=== Test 14: solve with trans=True ===")
b14 = ot.Point(20, 1.0)
x_normal = hodlr.solve(b14, False)
# The transposed solve is not implemented: trans=True must be rejected
with ott.assert_raises(RuntimeError):
    hodlr.solve(b14, True)
print("  PASS")

# === Test 15: 1D line, n=50 ===
print("\n=== Test 15: 1D line, n=50 ===")
vertices50 = ot.Sample([[i * 0.02] for i in range(50)])
af50 = HODLRTestAssemblyFunction(vertices50, 0.1)
params15 = make_params(leaf=25)
params15.setMaxRank(50)
hodlr50 = build_hodlr(vertices50, af50, params=params15)
K50 = build_dense_matrix(af50, 50)
b50 = ot.Point(50, 1.0)
ott.assert_almost_equal(hodlr50.solve(b50), K50.solveLinearSystem(b50), 1.0e-6, 1.0e-6)
print(f"  logDet= {hodlr50.logDeterminant():.6f}")
print("  PASS")

# === Test 16: gemv ===
print("\n=== Test 16: gemv ===")
x16 = [float(i) / 20 for i in range(20)]
x16_pt = ot.Point(x16)
y16 = ot.Point(20, 0.0)
hodlr.gemv('N', 1.0, x16_pt, 0.0, y16)
# Compare against dense matvec: y = K * x
y16_dense = ot.Point(20, 0.0)
for i in range(20):
    s = 0.0
    for j in range(20):
        s += af(i, j) * x16[j]
    y16_dense[i] = s
ott.assert_almost_equal(y16, y16_dense, 1.0e-6, 1.0e-6)
print("  PASS")

# === Test 17: gemv with alpha/beta ===
print("\n=== Test 17: gemv with alpha/beta ===")
x17 = [1.0] * 20
y17_init = [2.0] * 20
x17_pt = ot.Point(x17)
y17_pt = ot.Point(y17_init)
hodlr.gemv('N', 3.0, x17_pt, -1.0, y17_pt)
# y = 3*A*x - 1*y_init
y17_ref = ot.Point(20, 0.0)
for i in range(20):
    s = 0.0
    for j in range(20):
        s += af(i, j)
    y17_ref[i] = 3.0 * s - 2.0
ott.assert_almost_equal(y17_pt, y17_ref, 1.0e-10, 1.0e-10)
print("  PASS")

# === Test 18: HODLRMatrixFactory __repr__ ===
print("\n=== Test 18: HODLRMatrixFactory __repr__ ===")
factory = ot.HODLRMatrixFactory()
repr_f = repr(factory)
assert len(repr_f) > 0, "repr() should not be empty"
print("  PASS")

# === Test 19: HODLRMatrixParameters __repr__ / __str__ ===
print("\n=== Test 19: HODLRMatrixParameters __repr__ / __str__ ===")
assert len(repr(params)) > 0, "repr() should not be empty"
assert len(str(params)) > 0, "str() should not be empty"
print("  PASS")

# === Test 20: 2D grid, n=100 ===
print("\n=== Test 20: 2D grid, n=100 ===")
n20 = 9
intervalMesher = ot.IntervalMesher([n20, n20])
interval = ot.Interval([0.0] * 2, [1.0] * 2)
mesh2D = intervalMesher.build(interval)
vertices20 = mesh2D.getVertices()
size20 = vertices20.getSize()
af20 = HODLRTestAssemblyFunction(vertices20, 0.2)
params20 = make_params(leaf=25, factorization="LU")
params20.setMaxRank(50)
hodlr20 = build_hodlr(vertices20, af20, params=params20)
K20 = build_dense_matrix(af20, size20)
b20 = ot.Point(size20, 1.0)
ott.assert_almost_equal(hodlr20.solve(b20), K20.solveLinearSystem(b20), 1.0e-4, 1.0e-4)
print("  PASS")

# === Test 21: relative nugget via discretizeHODLRMatrix ===
print("\n=== Test 21: relative nugget ===")
nugget = 1.0e-4
ot.ResourceMap.SetAsScalar("HODLRMatrix-Nugget", nugget)
cov21 = ot.MaternModel([2.0], [1.0], 2.5)
grid21 = ot.Box([40]).generate()
grid21 *= 2.0
grid21 -= 1.0
n21 = grid21.getSize()
params21 = make_params(leaf=5, factorization="LLT")
params21.setMaxRank(40)
hodlr21 = cov21.discretizeHODLRMatrix(grid21, params21)
# kernel has unit amplitude: the diagonal grows by the relative nugget
diag_ref21 = ot.Point(n21, 1.0 + nugget)
ott.assert_almost_equal(hodlr21.getDiagonal(), diag_ref21, 1.0e-9, 1.0e-9)
hodlr21.factorize("LLT")
b21 = ot.Point(n21, 1.0)
K21 = ot.CovarianceMatrix(cov21.discretize(grid21))
for i in range(n21):
    K21[i, i] += nugget
ott.assert_almost_equal(hodlr21.solve(b21), K21.solveLinearSystem(b21), 1.0e-6, 1.0e-6)
ld21 = hodlr21.logDeterminant()
ld_ref21 = K21.computeLogAbsoluteDeterminant()[0]
ott.assert_almost_equal(ot.Point([ld21]), ot.Point([ld_ref21]), 1.0e-4, 1.0e-4)

# direct applyNugget on a factory-built matrix
ot.ResourceMap.SetAsScalar("HODLRMatrix-Nugget", 2.0e-3)
factory21 = ot.HODLRMatrixFactory()
hodlr21b = factory21.build(grid21, 1, True, make_params(leaf=5))
hodlr21b.assembleReal(HODLRTestAssemblyFunction(grid21, 0.1), 'L')
hodlr21b.applyNugget()
ott.assert_almost_equal(hodlr21b.getDiagonal(), ot.Point(n21, 1.0 + 2.0e-3), 1.0e-9, 1.0e-9)

# a zero nugget disables the regularization
ot.ResourceMap.SetAsScalar("HODLRMatrix-Nugget", 0.0)
hodlr21c = cov21.discretizeHODLRMatrix(grid21, params21)
ott.assert_almost_equal(hodlr21c.getDiagonal(), ot.Point(n21, 1.0), 1.0e-9, 1.0e-9)
ot.ResourceMap.SetAsScalar("HODLRMatrix-Nugget", 1.0e-8)
print("  PASS")

print("\n=== ALL TESTS PASSED ===")
