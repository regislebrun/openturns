/* Profiling driver for the hmat-oss backend: same workload as profile_ot. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int bench_hmat(const int n, const int dim, const double * points,
               const double length_scale, const double epsilon,
               const double nugget, const int leaf_size,
               const int factorization, const int compression,
               const double * b, double * x, double * logdet_out,
               double * t_assemble, double * t_factorize, double * t_solve,
               size_t * compressed, size_t * uncompressed, size_t * full_size);

int main(int argc, char ** argv)
{
  const int n = (argc > 1) ? atoi(argv[1]) : 1024;
  const int leaf = (argc > 2) ? atoi(argv[2]) : 64;
  const double eps = (argc > 3) ? atof(argv[3]) : 1.0e-6;
  const double l = (argc > 4) ? atof(argv[4]) : 0.1;

  const int m = (int)sqrt((double)n);
  double * pts = (double *) malloc((size_t) n * 2 * sizeof(double));
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < m; ++j)
    {
      pts[2 * (i * m + j) + 0] = i / (double)(m - 1);
      pts[2 * (i * m + j) + 1] = j / (double)(m - 1);
    }
  double * b = (double *) malloc((size_t) n * sizeof(double));
  double * x = (double *) malloc((size_t) n * sizeof(double));
  for (int i = 0; i < n; ++i) b[i] = 1.0;

  double logdet, t_asm, t_fact, t_solve;
  size_t comp, uncomp, full;
  int rc = bench_hmat(n, 2, pts, l, eps, 1.0e-6, leaf, 1 /*HODLRSYM*/,
                      3 /*aca_plus*/, b, x, &logdet,
                      &t_asm, &t_fact, &t_solve, &comp, &uncomp, &full);
  if (rc)
    fprintf(stderr, "bench_hmat failed rc=%d\n", rc);
  else
    printf("n=%d leaf=%d eps=%.1e: assemble=%.4fs factorize=%.4fs ratio=%.3f\n",
           n, leaf, eps, t_asm, t_fact, (double)comp / (double)full);
  free(pts); free(b); free(x);
  return rc;
}
