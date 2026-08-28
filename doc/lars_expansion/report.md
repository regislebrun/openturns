# LARS and OMP Expansion: Comparison Report

## Overview

This report compares the C++ `SparseExpansion` class (from PR #2987)
with the Python reference implementations of OMP and LARS by mbaudin47,
and with the Sobol index values reported in the PR discussion. It also
evaluates the impact of the `SparseExpansion-ConsecutiveIncreases`
stopping criterion.

## Test Case Parameters

Both test cases use the **Ishigami sparse problem**:

- **Model**: f(xi1, xi2, xi3) = sin(xi1) + 7\*(sin(xi2))^2 + 0.1\*xi3^4\*sin(xi1)
- **Distribution**: Uniform(-pi, pi)^3
- **Basis**: Legendre product polynomials with linear enumeration
- **Sample size**: n=75
- **Basis size**: 286 (max polynomial degree=10)
- **Reference Sobol S1**: X1=0.3139, X2=0.4424, X3=0.0000

The two test cases proposed by mbaudin47 (PR #2987, May 5, 2026) are:

1. **OMP**: Python implementation at
   https://gist.github.com/mbaudin47/052b9de72b5fa7e09d105151029a55a5
   (Algorithm B.1 from Luthen et al. 2021, page 628)
2. **LARS**: Python implementation at
   https://gist.github.com/mbaudin47/19cb6df17f1ef8cc3f18c50a4fb0e45e
   (Least Angle Regression Stepwise)

Both produce KFold cross-validation score plots as functions of iteration count.

## mbaudin47's Reported Sobol Index Values (PR #2987)

| Method | Fitting  | X1 S1   | X2 S1   | X3 S1 |
|--------|----------|---------|---------|-------|
| OMP    | KFold    | 0.3145  | 0.4415  | ~0    |
| OMP    | CLOO     | 0.3145  | 0.4415  | ~0    |
| LARS   | KFold    | 0.3144  | 0.4431  | ~0    |
| LARS   | CLOO     | 0.3149  | 0.4412  | ~0    |

Reference: X1 S1=0.314, X2 S1=0.442, X3 S1=0.0.

## Bug Fix: LARS Weighted Correlations

A bug was found in the LARS correlation computation. The code computed
`c_j = PhiW^T * r = Phi^T sqrt(W) r` instead of the correct weighted
correlation `c_j = <phi_j, r>_W = Phi^T W r`. The same issue affected
the `d` vector used in the step-size catch-up condition.

For uniform weights (`w_s = 1/n`), `sqrt(w_s) = 1/sqrt(n)` differs
from `w_s = 1/n` by a constant factor, so the signs were identical and
the bug was hidden. For non-uniform weights, the relative importance of
sample points changes, potentially flipping correlation signs and
selecting wrong basis functions.

The fix precomputes `sqrt(w_s)` and multiplies the residual and
prediction direction by it before correlating with `PhiW`, yielding the
correct `Phi^T W r` and `Phi^T W u`.

## Impact of Consecutive-Increases Stopping Criterion

The `SparseExpansion-ConsecutiveIncreases` ResourceMap key controls how
many consecutive increases in CV error are tolerated before stopping.
Results for k=1 (first increase), k=2 (default), and k=5:

| k | Method | Fitting | Active | X1 S1    | X2 S1    | X3 S1  | Max diff |
|---|--------|---------|--------|----------|----------|--------|----------|
| 1 | OMP    | KFold   | 26     | 0.31448  | 0.44151  | 0.0000 | 9.0e-04  |
| 1 | OMP    | CLOO    | 24     | 0.31447  | 0.44156  | 0.0000 | 8.6e-04  |
| 1 | LARS   | KFold   | 10     | 0.33469  | 0.43264  | 0.0004 | 2.1e-02  |
| 1 | LARS   | CLOO    | 7      | 0.41175  | 0.43694  | 0.0055 | 9.8e-02  |
| 2 | OMP    | KFold   | 34     | 0.31448  | 0.44151  | 0.0000 | 9.0e-04  |
| 2 | OMP    | CLOO    | 26     | 0.31448  | 0.44151  | 0.0000 | 9.0e-04  |
| 2 | LARS   | KFold   | 10     | 0.33469  | 0.43264  | 0.0004 | 2.1e-02  |
| 2 | LARS   | CLOO    | 10     | 0.33469  | 0.43264  | 0.0004 | 2.1e-02  |
| 5 | OMP    | KFold   | 41     | 0.31448  | 0.44153  | 0.0000 | 8.9e-04  |
| 5 | OMP    | CLOO    | 34     | 0.31448  | 0.44151  | 0.0000 | 9.0e-04  |
| 5 | LARS   | KFold   | 46     | 0.31483  | 0.44169  | 0.0000 | 9.2e-04  |
| 5 | LARS   | CLOO    | 10     | 0.33469  | 0.43264  | 0.0004 | 2.1e-02  |

### Key Observations

1. **OMP is robust to k**: All values of k produce excellent Sobol indices
   (max diff < 1e-3) because OMP's CV error decreases nearly monotonically.
   k=1 selects 26 functions, k=2 selects 34, k=5 selects 41.

2. **LARS with KFold and k=5 matches OMP**: With enough tolerance for CV
   oscillation, LARS KFold selects 46 functions and achieves
   max_diff=9.2e-04, comparable to OMP. This confirms the LARS algorithm
   is correct; it just needs more iterations to converge due to oscillating
   CV error.

3. **LARS with CLOO is limited**: CLOO produces higher CV errors than
   KFold, causing the `MaximumErrorFactor` stopping criterion to trigger
   early (at 10 active functions for k >= 2). This is a fitting-algorithm
   limitation, not an LARS algorithm issue.

4. **The default k=2 is a good compromise for OMP**: It allows one
   "bounce" in the CV error curve, which is sufficient for OMP.

## C++ SparseExpansion vs. Python Reference

### OMP Comparison (CLOO, k=2)

| Implementation       | Active | X1 S1    | X2 S1    | X3 S1  |
|----------------------|--------|----------|----------|--------|
| C++ OMP (k=2)        | 26     | 0.31448  | 0.44151  | 0.0000 |
| Python OMP (no stop) | 75     | 0.0203   | 0.0572   | 0.0438 |
| Reference            | --     | 0.3139   | 0.4424   | 0.0000 |

The Python OMP reference has **no CV stopping** and overfits badly when
run to completion (75 active = sample size).

### LARS Comparison (CLOO, k=2)

| Implementation        | Active | X1 S1    | X2 S1    | X3 S1  |
|-----------------------|--------|----------|----------|--------|
| C++ LARS (k=2)        | 10     | 0.33469  | 0.43264  | 0.0004 |
| Python LARS (no stop) | 75     | 0.0746   | 0.1374   | 0.0058 |
| Reference             | --     | 0.3139   | 0.4424   | 0.0000 |

## Differences Explained

### 1. CV Stopping vs. No Stopping

mbaudin47's Python implementations have **no cross-validation stopping**.
They run through all `basisSize - 1` iterations unconditionally, selecting
75 functions with 75 samples. This causes massive overfitting.

The C++ `SparseExpansion` uses CV-based early stopping with the
`SparseExpansion-ConsecutiveIncreases` key (default: 2), plus hard limits:

- `LeastSquaresMetaModelSelection-MaximumErrorFactor` (default 2.0)
- `LeastSquaresMetaModelSelection-MaximumError` (default 0.5)
- `LeastSquaresMetaModelSelection-ErrorThreshold` (default 0.0)

### 2. OMP vs. LARS Selection Path

- **OMP**: Picks the basis function with largest absolute correlation to the
  residual. Adds it and re-solves the LS system. The residual is fully
  updated after each step. CV error decreases nearly monotonically.
- **LARS**: Computes an equiangular direction and steps along it until
  another inactive function's correlation catches up. The "partial steps"
  cause the active set to grow differently, and CV error oscillates more.

### 3. Weighted Correlations

Both OMP and LARS use **weighted correlations** (`c = Phi^T * W * r`),
consistent with non-uniform quadrature weights. mbaudin47's Python
implementations use unweighted correlations (`c = X^T * r / n`).
With uniform weights (1/n), both are equivalent.

### 4. Why LARS KFold Needs k=5 to Match OMP

LARS equiangular steps cause the active set to grow in a different order
than OMP. The CV error oscillates (increases at some steps, decreases at
others), so k=2 stops LARS before it reaches the same accuracy as OMP.
With k=5, LARS tolerates more oscillation and eventually selects a
comparable active set (46 vs 41 for OMP).

## Summary

| Comparison                             | Result |
|----------------------------------------|--------|
| C++ OMP vs. mbaudin47 OMP values      | Match within 9e-4 (all k) |
| C++ LARS KFold (k=5) vs. mbaudin47    | Match within 9.2e-4 |
| C++ LARS KFold (k=2) vs. mbaudin47    | Diverges (2.1e-2) due to early CV stop |
| k=1 vs. k=5 for OMP                   | 26 vs 41 active functions, both accurate |
| k=1 vs. k=5 for LARS KFold            | 10 vs 46 active, k=5 matches OMP |
| Python ref vs. C++ (both methods)     | Python overfits (75 active), C++ avoids via CV stopping |

### Key Takeaways

1. The C++ `SparseExpansion` OMP faithfully reproduces mbaudin47's values.
2. LARS with KFold and k=5 also reproduces mbaudin47's values (within
   9.2e-4), confirming the algorithm is correct.
3. LARS with k=2 diverges because CV oscillation triggers early stopping;
   this is a stopping-criterion tuning issue, not an algorithm bug.
4. The consecutive-increases criterion (k=2 default) works well for OMP
   but is too aggressive for LARS with CLOO fitting.
5. The Python reference implementations lack CV stopping and overfit.

## Reproduction

All results are produced by the following Python tests (run from the
repository root with the installed library on PYTHONPATH):

```bash
# Ishigami sparse: Sobol indices for k=2 (default), KFold and CLOO
PYTHONPATH=$HOME/tmp/lars_expansion/install/lib64/python3.13/site-packages \
  python3 python/test/t_SparseExpansion_ishigami_sparse.py

# OMP comparison: C++ vs Python reference (CLOO, k=2)
PYTHONPATH=$HOME/tmp/lars_expansion/install/lib64/python3.13/site-packages \
  python3 python/test/t_SparseExpansion_omp_comparison.py

# LARS comparison: C++ vs Python reference (CLOO, k=2)
PYTHONPATH=$HOME/tmp/lars_expansion/install/lib64/python3.13/site-packages \
  python3 python/test/t_SparseExpansion_lars_comparison.py

# Consecutive-increases k sweep (k=1, k=2, k=5)
PYTHONPATH=$HOME/tmp/lars_expansion/install/lib64/python3.13/site-packages \
  python3 python/test/t_SparseExpansion_std.py
```

The C++ tests are registered in ctest:

```bash
cd $HOME/tmp/lars_expansion
ctest -R pyinstallcheck_SparseExpansion -V
```
