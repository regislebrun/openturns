// HODLR performance profiling driver for the OpenTURNS backend.
//
// Builds a symmetric exponential covariance matrix on a 2D grid, assembles it
// as a HODLR matrix and factorizes it, exactly like the benchmark driver does.
// It is used for callgrind profiling and for wall-clock timing without the
// Python interpreter.
//
// Compile with:
//   c++ -O2 -std=c++14 profile_ot.cxx -o profile_ot \
//       -I<repo>/build/install/include -L<repo>/build/install/lib64 -lOT \
//       -Wl,-rpath,<repo>/build/install/lib64
#include <openturns/OT.hxx>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>

using namespace OT;

static double now_s()
{
  timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + 1.0e-9 * ts.tv_nsec;
}

int main(int argc, char ** argv)
{
  const int n = (argc > 1) ? std::atoi(argv[1]) : 1600;
  const int leafSize = (argc > 2) ? std::atoi(argv[2]) : 64;
  const Scalar epsilon = (argc > 3) ? std::atof(argv[3]) : 1.0e-6;
  const Scalar lengthScale = (argc > 4) ? std::atof(argv[4]) : 0.1;
  const int iterations = (argc > 5) ? std::atoi(argv[5]) : 1;

  const int m = static_cast<int>(std::sqrt(static_cast<double>(n)));
  Sample vertices(0, 2);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < m; ++j)
      vertices.add(Point({i / Scalar(m - 1), j / Scalar(m - 1)}));

  MaternModel matern({lengthScale}, {1.0}, 0.5);
  IsotropicCovarianceModel model(matern, 2);

  HODLRMatrixParameters params;
  params.setAssemblyEpsilon(epsilon);
  params.setRecompressionEpsilon(epsilon);
  params.setMinLeafSize(leafSize);
  params.setFactorizationMethod("LLT");

  ResourceMap::SetAsScalar("HODLRMatrix-Nugget", 1.0e-6);

  double tAssemble = 0.0;
  double tFactorize = 0.0;
  for (int it = 0; it < iterations; ++it)
  {
    double t0 = now_s();
    HODLRMatrix hodlr = model.discretizeHODLRMatrix(vertices, params);
    double t1 = now_s();
    hodlr.factorize("LLT");
    double t2 = now_s();
    if (it > 0 || iterations == 1)
    {
      tAssemble += t1 - t0;
      tFactorize += t2 - t1;
    }
  }
  std::printf("n=%d leaf=%d eps=%.1e: assemble=%.4fs factorize=%.4fs (avg over %d runs)\n",
              n, leafSize, epsilon, tAssemble, tFactorize, iterations);
  return 0;
}
