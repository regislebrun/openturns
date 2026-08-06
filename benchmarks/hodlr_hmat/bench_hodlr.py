#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Benchmark comparing the OpenTURNS HODLR backend (Base/Stat/HODLR*) with the
HODLR backend of the hmat-oss library (https://github.com/jeromerobert/hmat-oss).

Both backends assemble the same symmetric covariance matrix
    K(i, j) = exp(-|p_i - p_j| / length_scale) + nugget * delta_{i,j}
on the same set of points, factorize it and solve the same linear system
A x = b. The benchmark reports wall-clock times (assembly, factorization,
solve), the compression ratio, the log-determinant and the accuracy of the
solution and of the log-determinant w.r.t. a dense reference (when the dense
matrix fits in memory).

OpenTURNS side:  ot.MaternModel([l], [1.0], 0.5).discretizeHODLRMatrix()
                 (Matern nu=0.5 is exactly the exponential kernel). Assembly
                 runs in native C++ through HODLRCovarianceAssemblyFunction,
                 the same path as GaussianProcess.
hmat-oss side:   libbench_hmat.so, a thin C shim (bench_hmat.c) built on the
                 hmat-oss C API.

Usage:
  python3 bench_hodlr.py [--sizes 1000,2500,5000]
                         [--epsilons 1e-4,1e-6,1e-8]
                         [--leafs 16,64]
                         [--points grid2d|grid1d|random2d]
                         [--length 0.1]
                         [--nugget 1e-6]
                         [--method LLT|LU]
                         [--compression aca_plus|aca_partial|aca_full|aca_random|svd|rrqr]
                         [--threads 1]
                         [--dense-limit 1073741824]
                         [--out results.json]
                         [--no-rebuild]
"""
import argparse
import ctypes
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

# Control BLAS/OpenMP threading *before* importing numpy/OpenTURNS.
THREADS = os.environ.get("OPENBLAS_NUM_THREADS", "1")
for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[var] = str(os.environ.get(var, "1"))
os.environ["OMP_DYNAMIC"] = "FALSE"

import openturns as ot  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SHIM_SRC = os.path.join(HERE, "bench_hmat.c")
SHIM_OUT = os.path.join(REPO_ROOT, "build", "benchmarks", "libbench_hmat.so")

COMPRESSIONS = {
    "svd": 0,
    "aca_full": 1,
    "aca_partial": 2,
    "aca_plus": 3,
    "aca_random": 4,
    "rrqr": 5,
}

# Map an OpenTURNS factorization method to the hmat-oss one.
# hmat-oss HODLRSYM logdet returns log|det(L)|, i.e. half of log det A; the
# driver doubles it so both backends report log det A.
HMAT_FACT = {"LLT": ("hodlrsym", 1, 2.0), "LU": ("hodlr", 0, 1.0)}


def build_shim(force=False):
    if not force and os.path.exists(SHIM_OUT):
        return SHIM_OUT
    os.makedirs(os.path.dirname(SHIM_OUT), exist_ok=True)
    cmd = [
        "cc", "-O3", "-fPIC", "-shared", SHIM_SRC, "-o", SHIM_OUT,
        "-I/usr/local/include", "-L/usr/local/lib", "-Wl,-rpath,/usr/local/lib",
        "-lhmat", "-lm",
    ]
    subprocess.check_call(cmd)
    return SHIM_OUT


class HmatBench:
    """ctypes binding of the bench_hmat shim."""

    def __init__(self, shim):
        self.lib = ctypes.CDLL(shim)
        f = self.lib.bench_hmat
        f.argtypes = [
            ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_double,
            ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int,
            ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
            ctypes.c_void_p, ctypes.c_void_p,
        ]
        f.restype = ctypes.c_int
        self.f = f

    def run(self, points, length, eps, nugget, leaf, hmat_fact,
            compression, b):
        n = len(points)
        x = np.zeros(n)
        logdet = ctypes.c_double(0.0)
        t_assemble = ctypes.c_double(0.0)
        t_factorize = ctypes.c_double(0.0)
        t_solve = ctypes.c_double(0.0)
        compressed = ctypes.c_size_t(0)
        uncompressed = ctypes.c_size_t(0)
        full = ctypes.c_size_t(0)
        rc = self.f(
            n, points.shape[1],
            points.ctypes.data_as(ctypes.c_void_p),
            length, eps, nugget, leaf,
            hmat_fact[1], compression,
            b.ctypes.data_as(ctypes.c_void_p),
            x.ctypes.data_as(ctypes.c_void_p),
            ctypes.byref(logdet),
            ctypes.byref(t_assemble), ctypes.byref(t_factorize),
            ctypes.byref(t_solve),
            ctypes.byref(compressed), ctypes.byref(uncompressed),
            ctypes.byref(full),
        )
        if rc:
            raise RuntimeError(f"bench_hmat returned {rc}")
        return {
            "t_assemble": t_assemble.value,
            "t_factorize": t_factorize.value,
            "t_solve": t_solve.value,
            "logdet": hmat_fact[2] * logdet.value,
            "x": x,
            "compressed": compressed.value,
            "uncompressed": uncompressed.value,
            "full_size": full.value,
            "compression_ratio": compressed.value / max(1, uncompressed.value),
        }


def make_points(kind, n):
    """Return (points (n x d), actual_n) for the requested geometry."""
    rng = np.random.default_rng(0)
    if kind == "grid1d":
        pts = np.linspace(0.0, 1.0, n)[:, None]
    elif kind == "grid2d":
        m = max(2, int(math.floor(math.sqrt(n))))
        xx, yy = np.meshgrid(np.linspace(0.0, 1.0, m), np.linspace(0.0, 1.0, m))
        pts = np.stack([xx.ravel(), yy.ravel()], axis=1)
    elif kind == "random2d":
        pts = rng.random((n, 2))
    else:
        raise ValueError(f"unknown --points {kind}")
    return np.ascontiguousarray(pts, dtype=np.float64)


def make_rhs(points):
    """Deterministic smooth right-hand side."""
    s = points.sum(axis=1)
    return np.ascontiguousarray(np.sin(np.pi * s), dtype=np.float64)


def dense_reference(points, length, nugget, b):
    n = len(points)
    diff = points[:, None, :] - points[None, :, :]
    d = np.linalg.norm(diff, axis=2)
    k = np.exp(-d / length) + nugget * np.eye(n)
    x = np.linalg.solve(k, b)
    sign, logdet = np.linalg.slogdet(k)
    return x, logdet


def run_ot(points, length, eps, nugget, leaf, method):
    ot.ResourceMap.SetAsScalar("HODLRMatrix-Nugget", nugget)
    params = ot.HODLRMatrixParameters()
    params.setAssemblyEpsilon(eps)
    params.setRecompressionEpsilon(eps)
    params.setMinLeafSize(leaf)
    params.setFactorizationMethod(method)
    # Matern nu=0.5 is the exponential kernel; the isotropic wrapper gives the
    # model the input dimension of the points (MaternModel scale is 1D).
    matern = ot.MaternModel([length], [1.0], 0.5)
    model = ot.IsotropicCovarianceModel(matern, points.shape[1])
    b = ot.Point(make_rhs(points))

    t0 = time.perf_counter()
    hodlr = model.discretizeHODLRMatrix(ot.Sample(points), params)
    t1 = time.perf_counter()
    hodlr.factorize(method)
    t2 = time.perf_counter()
    x = np.array(hodlr.solve(b))
    t3 = time.perf_counter()
    logdet = hodlr.logDeterminant()
    compressed, total = hodlr.compressionRatio()

    return {
        "t_assemble": t1 - t0,
        "t_factorize": t2 - t1,
        "t_solve": t3 - t2,
        "logdet": logdet,
        "x": x,
        "compression_ratio": compressed / float(max(1, total)),
        "full_size": total,
    }


def rel_err(a, b):
    denom = np.linalg.norm(b)
    return float(np.linalg.norm(a - b) / denom) if denom else 0.0


def jsonable(obj):
    """Make an object JSON-serializable (numpy arrays -> lists)."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    return obj


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sizes", default="625,1600,2500",
                    help="matrix sizes, comma separated")
    ap.add_argument("--epsilons", default="1e-4,1e-6,1e-8",
                    help="assembly/compression epsilons")
    ap.add_argument("--leafs", default="16,64",
                    help="minimum leaf sizes")
    ap.add_argument("--points", default="grid2d", choices=["grid1d", "grid2d", "random2d"])
    ap.add_argument("--length", type=float, default=0.1,
                    help="kernel length scale (correlation length)")
    ap.add_argument("--nugget", type=float, default=1.0e-6,
                    help="diagonal nugget added on both sides")
    ap.add_argument("--method", default="LLT", choices=["LLT", "LU"])
    ap.add_argument("--compression", default="aca_plus",
                    choices=sorted(COMPRESSIONS))
    ap.add_argument("--threads", type=int, default=1,
                    help="BLAS/OpenMP threads for both backends")
    ap.add_argument("--dense-limit", type=float, default=1.0e9,
                    help="build the dense reference only if n^2*8 bytes is below this")
    ap.add_argument("--out", default="", help="write a JSON report to this file")
    ap.add_argument("--no-rebuild", action="store_true",
                    help="do not rebuild the C shim if it exists")
    args = ap.parse_args(argv)

    for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = str(args.threads)

    sizes = [int(s) for s in args.sizes.split(",")]
    epsilons = [float(e) for e in args.epsilons.split(",")]
    leafs = [int(v) for v in args.leafs.split(",")]
    hmat_fact = HMAT_FACT[args.method]
    compression = COMPRESSIONS[args.compression]

    shim = build_shim(force=not args.no_rebuild)
    hmat = HmatBench(shim)

    results = []
    print("n       eps       leaf  backend  t_asm[s]  t_fact[s]  t_solve[s]  ratio   "
          "err_solve  err_logdet")
    for n in sizes:
        points = make_points(args.points, n)
        b = make_rhs(points)
        dense_x = None
        dense_logdet = None
        if n * n * 8.0 < args.dense_limit:
            dense_x, dense_logdet = dense_reference(points, args.length, args.nugget, b)
        for eps in epsilons:
            for leaf in leafs:
                ot_res = run_ot(points, args.length, eps, args.nugget, leaf,
                                args.method)
                hmat_res = hmat.run(points, args.length, eps, args.nugget, leaf,
                                    hmat_fact, compression, b)

                entries = [("ot", ot_res), ("hmat", hmat_res)]
                if dense_x is not None:
                    per_backend = [
                        (tag, rel_err(res["x"], dense_x),
                         abs(res["logdet"] - dense_logdet) / abs(dense_logdet))
                        for tag, res in entries
                    ]
                else:
                    cross = rel_err(ot_res["x"], hmat_res["x"])
                    ld_cross = abs(ot_res["logdet"] - hmat_res["logdet"]) / abs(hmat_res["logdet"])
                    per_backend = [("ot", cross, ld_cross), ("hmat", cross, ld_cross)]

                for tag, res, err, ld_err in (
                    (t, r, e, l) for (t, r), (e, l) in
                    zip(entries, [(e, l) for _, e, l in per_backend])
                ):
                    print(f"{n:6d} {eps:8.1e} {leaf:6d}  {tag:<7s} "
                          f"{res['t_assemble']:8.4f} {res['t_factorize']:8.4f} "
                          f"{res['t_solve']:9.5f}  {res['compression_ratio']:6.3f} "
                          f"{err:9.2e} {ld_err:10.2e}")

                results.append({
                    "n": n, "eps": eps, "leaf": leaf, "points": args.points,
                    "length": args.length, "nugget": args.nugget,
                    "method": args.method, "compression": args.compression,
                    "ot": ot_res, "hmat": hmat_res,
                    "err_solve": per_backend[0][1],
                    "err_logdet": per_backend[0][2],
                })

    if args.out:
        with open(args.out, "w") as f:
            json.dump(jsonable(results), f, indent=2)
        print(f"\nreport written to {args.out}")


if __name__ == "__main__":
    sys.exit(main())
