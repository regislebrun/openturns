/*
  HODLR benchmark: thin C wrapper around the hmat-oss C API.

  It assembles the symmetric covariance matrix K(i, j) = exp(-|p_i - p_j| / l)
  on a set of points, factorizes it as a HODLR matrix and solves one
  linear system, exposing timings and compression statistics to the
  Python driver (bench_hodlr.py) through ctypes.

  Compile with:
    cc -O3 -fPIC -shared bench_hmat.c -o libbench_hmat.so \
       -I/usr/local/include -L/usr/local/lib -lhmat -lm
*/
#define _GNU_SOURCE
#include <hmat/hmat.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct
{
  int n;
  int dim;
  const double * points;
  double length_scale;
  double nugget;
} bench_data_t;

static double bench_now()
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + 1.0e-9 * ts.tv_nsec;
}

static void simple_compute(void * vdata, int i, int j, void * result)
{
  const bench_data_t * d = (const bench_data_t *) vdata;
  const double * pi = d->points + (size_t) i * d->dim;
  const double * pj = d->points + (size_t) j * d->dim;
  double acc = 0.0;
  for (int k = 0; k < d->dim; ++k)
  {
    const double dx = pi[k] - pj[k];
    acc += dx * dx;
  }
  *((double *) result) = exp(-sqrt(acc) / d->length_scale) + (i == j ? d->nugget : 0.0);
}

/* factorization: 0 -> HODLR (LU), 1 -> HODLRSYM (symmetric Cholesky-like)
   compression: index into hmat_compress_t
     (0 svd, 1 aca_full, 2 aca_partial, 3 aca_plus, 4 aca_random, 5 rrqr)
   b: right-hand side in original numbering (size n)
   x: solution in original numbering (size n, output)
   All timings are in seconds, all sizes are counts of stored coefficients. */
int bench_hmat(const int n,
               const int dim,
               const double * points,
               const double length_scale,
               const double epsilon,
               const double nugget,
               const int leaf_size,
               const int factorization,
               const int compression,
               const double * b,
               double * x,
               double * logdet_out,
               double * t_assemble,
               double * t_factorize,
               double * t_solve,
               size_t * compressed,
               size_t * uncompressed,
               size_t * full_size)
{
  hmat_interface_t hmat;
  hmat_init_default_interface(&hmat, HMAT_DOUBLE_PRECISION);

  bench_data_t data = {n, dim, points, length_scale, nugget};

  hmat_clustering_algorithm_t * ctm = hmat_create_clustering_median();
  hmat_clustering_algorithm_t * ctmd = hmat_create_clustering_max_dof(ctm, leaf_size);
  hmat_cluster_tree_t * tree = hmat_create_cluster_tree((double *) points, dim, n, ctmd);
  hmat_delete_clustering(ctmd);
  hmat_delete_clustering(ctm);
  if (!tree)
  {
    fprintf(stderr, "bench_hmat: cluster tree creation failed\n");
    return 1;
  }

  hmat_admissibility_t * ac = hmat_create_admissibility_hodlr();
  hmat_matrix_t * matrix = hmat.create_empty_hmatrix_admissibility(tree, tree, 1, ac);
  hmat_delete_admissibility(ac);
  if (!matrix)
  {
    fprintf(stderr, "bench_hmat: empty hmatrix creation failed\n");
    return 2;
  }
  hmat.own_cluster_trees(matrix, 1, 0);
  hmat.set_low_rank_epsilon(matrix, epsilon);

  memcpy(x, b, (size_t) n * sizeof(double));
  hmat.vector_reorder(x, tree, 0, NULL, 1);

  hmat_compression_algorithm_t * comp = NULL;
  switch (compression)
  {
    case 0: comp = hmat_create_compression_svd(epsilon); break;
    case 1: comp = hmat_create_compression_aca_full(epsilon); break;
    case 2: comp = hmat_create_compression_aca_partial(epsilon); break;
    case 3: comp = hmat_create_compression_aca_plus(epsilon); break;
    case 4: comp = hmat_create_compression_aca_random(epsilon); break;
    case 5: comp = hmat_create_compression_rrqr(epsilon); break;
    default: comp = hmat_create_compression_aca_plus(epsilon); break;
  }
  if (!comp)
  {
    fprintf(stderr, "bench_hmat: compression creation failed\n");
    return 6;
  }

  hmat_assemble_context_t acx;
  hmat_assemble_context_init(&acx);
  acx.compression = comp;
  acx.user_context = &data;
  acx.simple_compute = simple_compute;
  acx.lower_symmetric = 1;

  double t0 = bench_now();
  const int rc_assemble = hmat.assemble_generic(matrix, &acx);
  double t1 = bench_now();
  hmat_delete_compression(comp);
  if (rc_assemble)
  {
    fprintf(stderr, "bench_hmat: assemble_generic failed (%d)\n", rc_assemble);
    return 3;
  }
  *t_assemble = t1 - t0;

  hmat_factorization_context_t fcx;
  hmat_factorization_context_init(&fcx);
  fcx.factorization = factorization ? hmat_factorization_hodlrsym : hmat_factorization_hodlr;
  t0 = bench_now();
  const int rc_factorize = hmat.factorize_generic(matrix, &fcx);
  t1 = bench_now();
  if (rc_factorize)
  {
    fprintf(stderr, "bench_hmat: factorize_generic failed (%d)\n", rc_factorize);
    return 4;
  }
  *t_factorize = t1 - t0;

  t0 = bench_now();
  const int rc_solve = hmat.solve_dense(matrix, x, 1);
  t1 = bench_now();
  if (rc_solve)
  {
    fprintf(stderr, "bench_hmat: solve_dense failed (%d)\n", rc_solve);
    return 5;
  }
  hmat.vector_restore(x, tree, 0, NULL, 1);
  *t_solve = t1 - t0;

  double logdet = 0.0;
  if (hmat.logdet(matrix, &logdet))
    logdet = 0.0 / 0.0;
  *logdet_out = logdet;

  hmat_info_t info;
  memset(&info, 0, sizeof(info));
  hmat.get_info(matrix, &info);
  *compressed = info.compressed_size;
  *uncompressed = info.uncompressed_size;
  *full_size = info.full_size;

  hmat.destroy(matrix);
  return 0;
}
