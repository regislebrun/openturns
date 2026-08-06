# HODLR benchmark: OpenTURNS vs hmat-oss

This benchmark compares the HODLR backend of OpenTURNS
(`lib/src/Base/Stat/HODLR*`, in this `hodlr` branch) with the HODLR backend of
the [hmat-oss](https://github.com/jeromerobert/hmat-oss) library (v1.11.1 used
here).

Both backends assemble the **same** symmetric covariance matrix

    K(i, j) = exp(-|p_i - p_j| / length_scale) + nugget * delta(i, j)

on the **same** point set, factorize it (Cholesky / LLt) and solve the **same**
linear system A x = b. Reported metrics:

* wall-clock time of assembly, factorization and solve,
* compression ratio (stored coefficients / n^2),
* log-determinant,
* relative error of the solution and of the log-determinant with respect to a
  dense reference (computed with numpy when `n^2 * 8 bytes < --dense-limit`,
  otherwise the two backends are checked against each other).

## Files

| File | Role |
|------|------|
| `bench_hodlr.py` | Python driver (OpenTURNS via SWIG, hmat-oss via ctypes) |
| `bench_hmat.c` | Thin C shim around the hmat-oss C API (compiled on the fly) |

## Prerequisites

* a build of this OpenTURNS branch, with the Python bindings installed
  (`PYTHONPATH` must point at the site-packages directory),
* the hmat-oss library installed with its C header, e.g.
  `/usr/local/lib/libhmat.so` and `/usr/local/include/hmat/hmat.h`,
* a C compiler and numpy.

The C shim is compiled on first use into
`<repo>/build/benchmarks/libbench_hmat.so`; pass `--no-rebuild` to reuse it.

## Usage

```sh
PYTHONPATH=<repo>/build/install/lib64/python3.13/site-packages \
    python3 bench_hodlr.py --sizes 625,1600,2500 \
                           --epsilons 1e-4,1e-6,1e-8 \
                           --leafs 16,64 \
                           --points grid2d \
                           --method LLT \
                           --out results.json
```

Options:

* `--points grid2d|grid1d|random2d` point set geometry (domain `[0,1]^d`).
* `--length` kernel length scale (correlation length), default `0.1`.
* `--nugget` diagonal nugget added on both sides, default `1e-6`. A small
  nugget is needed because the pure exponential covariance is
  ill-conditioned, which makes the dense log-determinant meaningless; a nugget
  is also what OpenTURNS applies by default (`HODLRMatrix-Nugget`).
* `--method LLT|LU` factorization. `LLT` maps to `factorize("LLT")` on the
  OpenTURNS side and `hmat_factorization_hodlrsym` on the hmat-oss side; `LU`
  maps to `factorize("LU")` and `hmat_factorization_hodlr`.
* `--compression aca_plus|aca_partial|aca_full|aca_random|svd|rrqr` hmat-oss
  assembly compression (default `aca_plus`; `aca_random` is what hmat-oss's
  own `hodlrvsllt` example uses).
* `--threads` BLAS/OpenMP threads for both backends (default `1`; the
  OpenTURNS backend pins OpenBLAS to one thread anyway).
* `--out file.json` write the full report.

## Notes and caveats

* **Leaf size semantics differ.** OpenTURNS `HODLRMatrix-MinLeafSize` is a
  *minimum* leaf size (leaves end up in `[leaf, 2*leaf)`), while hmat-oss
  `hmat_create_clustering_max_dof` caps the leaf at the given size. The two
  trees are therefore not identical.
* **Ordering differs.** OpenTURNS reorders the points along a KDTree
  space-filling curve; hmat-oss uses a median clustering. Both are part of the
  measured pipeline.
* **log-determinant convention differs.** hmat-oss `HODLRSYM` factorization
  `logdet` returns `log|det(L)|`, i.e. **half** of `log det A` (its leaves are
  stored as Cholesky factors, see `src/hodlr.cpp` `logdet()`); the driver
  doubles it so both backends report `log det A`. For the `HODLR` (LU)
  factorization the value is already `log det A`.
* **Assembly epsilon maps to two knobs in hmat-oss:** `set_low_rank_epsilon`
  (recompression) and the compression algorithm epsilon (assembly). Both are
  set to the same value, mirroring OpenTURNS where
  `HODLRMatrix-AssemblyEpsilon` and `-RecompressionEpsilon` are aligned.
* **OpenTURNS assembly uses the native `HODLRCovarianceAssemblyFunction`**
  (the `discretizeHODLRMatrix` path, also used by `GaussianProcess`), not a
  per-entry Python callback, so the comparison is kernel-for-kernel.
* **Observed scaling.** On a `grid2d` set with `length=0.1`, the OpenTURNS
  factorization is currently much slower than hmat-oss at `n >= 1600`
  (seconds vs tens of milliseconds at `n=2500`); both track the requested
  epsilon in accuracy. This is an experimental backend: expect the numbers to
  move as the implementation matures.
