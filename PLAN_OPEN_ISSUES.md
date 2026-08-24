# Plan: Address the open issues and enhancements collected from the HODLR benchmark reports

## Context

Collected from the LaTeX reports in `benchmarks/hodlr_hmat/report/` and
`report_2026_08_11/` (5 reports), cross-checked against the current branch
state (HEAD `4b1d2bd99`, branch `hodlr`). Items already implemented
(nugget via `applyNugget`, adaptive rank, dense fallback for full-rank blocks,
space-filling ordering, recompression, log-det fix, P0-P3 perf fixes, P4
conditioning warning + `getRegularizationShift`) are NOT listed here except
where a status line below needs updating.

Priority order chosen by the user: performance first, then precision/robustness,
then benchmark metrics and MLE treatment.

## Open items (ranked)

### P1. HODLR LU factorization is 6-12x slower than its own Cholesky

**Source:** `benchmarks/hodlr_hmat/report/report.tex`, Findings paragraph "LU factorization":
native HODLR LU is 6-12x slower than its own Cholesky (223 ms at n=5000/1D, 705 ms at 3D),
while hmat-oss HODLR (LU) and the HMAT LU stay close to their symmetric cost.

**Status:** CLOSED. The HODLR backend now offers a single LLt (Cholesky)
factorization: the LU factorization path was removed from
`HODLRMatrixImplementation` / `HODLRCore`, `HODLRMatrixParameters` no longer has
a `FactorizationMethod` attribute, and `HODLRMatrix::factorize()` takes no
method argument. `KarhunenLoeveP1Algorithm`, GaussianProcess and the benchmark
drivers were switched to the Cholesky path, which also matches what
hmat-oss's HODLRSYM exposes.

### P2. Recompression stack dominates the factorization

**Source:** `report_2026_08_11/benchmark_sweep/report_sweep.md`, "Hot spots":
the low-rank recompression stack (QR + Gram-matrix DGEMM + syevd + back-projection)
accounts for ~77-82% of the HODLR factorization time and is "the next most profitable target".

**Status:** LARGELY CLOSED. Two commits cut the dominant recompression cost:
- `a7f19d3fe` (HODLRCore: apply the QR reflectors directly with dormqr): the
  back-projection phase went from O(nRows * k1^2) to O(nRows * k1 * rankOut);
  recompression 277 -> 231 ms on a 2D n=10000 Matern case.
- `4b1d2bd99`+ (HODLRCore: pre-truncate the correction stack by column weight,
  keys `HODLRMatrix-StackTruncation` / `-StackTruncationFactor`): the stacked
  corrections from all descendant levels are dropped before the QR+Gram+SVD when
  their cumulative weight `||u_j|| * ||v_j||` stays under `factor * tol *
  weightSum` (a provable bound via the triangle inequality, so the exact SVD on
  the reduced stack still meets the assembly tolerance). The factorization on
  the 2D 99x99 Matern benchmark went from ~1440 ms to ~670 ms (2.1x) with
  unchanged accuracy (solve/log-det match the dense reference to 1e-7/1e-10).

**Remaining (low priority):**
1. Skip the SVD recompression when the sum of correction ranks is small
   (totalRank <= rankIn + 1): the Ucat/Vcat concatenation is already nearly
   low-rank, truncate by removing the smallest correction instead.
2. Avoid copying U/V blocks into Ucat/Vcat (`HODLRCore.cxx:604-619`) by
   assembling the QR input via a strided view / block layout.
3. Partial eigendecomposition of the Gram (`syevr`/`dstegr`/`dsyevx` range='I')
   was evaluated and abandoned: OpenBLAS 0.3.30 ignores the range restriction
   (no speedup vs `dsyevd`) and its `dstegr`/`dstemr` corrupt the heap at
   n >= 200. Not pursued.
4. Verify the 2D/3D n=10000 factorization lands within ~1.5x of hmat-oss.

### P3. Solve gap vs HMAT (partially closed)

**Source:** `report_2026_08_11/rapport_hodlr_precision/` sec. "Performance":
the HODLR solve was 5-10x slower than HMAT (dominant cost `applyInverse`/substitution);
the sweep now reports 1.4-2.4x, so the gap is largely closed but remains.

**Status:** PARTIALLY CLOSED. Remaining gap is 1.4-2.4x in the sweep.

**Proposed work (low priority):**
1. Profile `applyInverse` / the solve substitution for a 2D n=20000 case.
2. Reduce per-node overhead: the reports point at `vector::operator[]`
   (~6.8%) and `vector::size()` (~6.7%) in the substitution hot loop;
   hoist size lookups and use raw pointers / contiguous buffers.
3. Only pursue if P1/P2 do not already absorb the gap.

### P4. Conditioning-aware failure signal

**Source:** `report_2026_08_11/hodlr_cholesky_accuracy/` recommendation 4 and
`hmat_1d_analysis/REPORT_HMAT_1D_ANALYSIS.tex` sec. 5.
The regularization loops (`HODLRMatrixImplementation.cxx:291`,
`HMatrixImplementation::factorize`) silently factorize `A + lambda I` when the
factorization fails; when lambda becomes non-negligible relative to the diagonal,
all post-factorization results (solve, log-determinant) are for a matrix far
from A, with no warning.

**Status:** CLOSED. Commit `13bbcbd82` makes both backends record the shift
actually applied (`getRegularizationShift()`) and emit a warning when it is
non-negligible relative to the diagonal, controlled by the
`HODLRMatrix-RegularizationWarnThreshold` / `HMatrix-RegularizationWarnThreshold`
resource keys. Covered by `t_HODLRMatrix_std.py` and `t_HMatrix_std.py`.

### P5. Accuracy metrics in the benchmark

**Source:** `report_2026_08_11/hodlr_cholesky_accuracy/` recommendation 5
(citing Litvinenko et al. 2019): report the KL divergence between the exact and
HODLR-approximated Gaussian laws and/or `||C C_H^-1 - I||`, so compression
quality is visible alongside runtime.

**Status:** OPEN.

**Proposed work:**
1. Extend `python/test/benchmark_gp_full.py` (and/or the sweep harness) with a
   `--accuracy` mode that, for n <= a dense limit, computes:
   - `||C - C_H||_2 / ||C||_2`
   - `||C C_H^-1 - I||_2`
   - `KL` (Litvinenko 2019 formula, via logdet and trace of `C C_H^-1`)
   using `discretizeHODLRMatrix` + the dense reference from LAPACK.
2. Add these as columns in the printed tables and/or the LaTeX sweep tables.
3. Document the metric definitions in the benchmark README.

**Status:** CLOSED. `python/test/benchmark_gp_full.py` gained an `--accuracy`
mode (`--accuracy-max-n` dense limit, default 2000; `--eps` assembly epsilon,
default 1e-6) that reports `||C - C_H||_2 / ||C||_2` (power iteration via HODLR
`gemv`), `||C C_H^-1 - I||_2` (largest singular value of `X - I`, with
`X = C C_H^-1` built by HODLR column solves) and the Litvinenko 2019 KL formula
via dense logdet, HODLR `logDeterminant()` and `tr(X)`, plus `log|C_H| - log|C|`
separately. The metric definitions are documented in
`benchmarks/hodlr_hmat/README.md`. The mode exposes the expected behaviour: at
`corr=0.1` the compression stays accurate to ~1e-8 at n~1000, while at
`corr=1.0` the near-singular regime drives the metrics (KL, logdet diff) to
large values at n >= 513, i.e. the regularized-factorization problem that P4
warns about.

### P6. Geoga et al. MLE treatment

**Source:** `report_2026_08_11/hodlr_cholesky_accuracy/` recommendation 6
(citing Geoga et al. 2020): for long length scales (theta ~ domain size) the
approximated likelihood surface is non-convex / loses positive definiteness;
the reference work uses a first-order optimizer, or one can tighten the
compression tolerance inside the fitter.

**Status:** OPEN.

**Proposed work:**
1. Investigate `GaussianProcessFitter` behaviour at corr=1.0 (already covered in
   `benchmark_gp_full.py`): check whether the HODLR backend currently converges
   to the LAPACK optimum.
2. If needed, expose a ResourceMap key
   (e.g. `GaussianProcessFitter-HODLRCompressionEpsilon`) letting the fitter
   tighten `HODLRMatrixParameters.assemblyEpsilon` during optimization.
3. Document the limitation in the GaussianProcessFitter docstring (Notes section)
   and in the ChangeLog.

## Implementation order

1. **P1 + P2** (both live in `recompressLowRank`/`compute`/`computeCholesky`
   in `HODLRCore.cxx`; profile first, then optimize). P1 closed, P2 largely
   closed by `a7f19d3fe` + the stack pre-truncation; only the low-priority
   sub-items above remain.
2. **P4** (conditioning-aware warning + `getRegularizationShift()` + test). CLOSED.
3. **P5** (benchmark accuracy metrics).
4. **P3** (solve gap) and **P6** (MLE treatment), low priority.

## Verification

- C++ build: `cmake --build build --target install --parallel $(( $(nproc) / 2 ))`
- Tests: `ctest -R "HODLRMatrix|GaussianProcess|HMatrix" -V`
- Perf regression: re-run `benchmarks/hodlr_hmat/bench_hodlr_vs_hmat.py` sweep
  and regenerate `report_sweep` tables; target factorize/solve within ~1.5-2x
  of hmat-oss.
- Accuracy must not regress: matvec relative error stays at the assembly epsilon
  on the 1D-4D sweep cases.
