//                                               -*- C++ -*-
/**
 *  @brief The test file of class GaussKronrod
 *
 *  Copyright 2005-2026 Airbus-EDF-IMACS-ONERA-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "openturns/OT.hxx"
#include "openturns/OTtestcode.hxx"
#include <cmath>

using namespace OT;
using namespace OT::Test;

typedef Collection<Complex> ComplexCollection;

int main(int, char *[])
{
  TESTPREAMBLE;
  OStream fullprint(std::cout);
  PlatformInfo::SetNumericalPrecision(8);

  try
  {
    const Point ref = {0, 0, 1.0, 0.3591409, 0.09390606, 0.01935569};
    for (UnsignedInteger n = 2; n < 6; ++ n)
    {
      String sum = "x0";
      for (UnsignedInteger i = 1; i < n; ++ i)
        sum += (OSS() << "+x" << i);
      SymbolicFunction f(Description::BuildDefault(n, "x"), Description({"exp(" + sum + ")"}));
      fullprint << "f=" << f.__str__() << std::endl;
      Sample vertices(1, n);
      for (UnsignedInteger j = 0; j < n; ++ j)
      {
        Point e(n);
        e[j] = 1.0;
        vertices.add(e);
      }
      Indices indices(n + 1);
      indices.fill();
      IndicesCollection simplicies(Collection<Indices>(1, indices));
      Mesh canonicalSimplex(vertices, simplicies, true);
      for (UnsignedInteger k = 1; k < 5; ++ k)
      {
        SimplicialCubature algo;
        algo.setRule(k);
        // fullprint << "Algo=" << algo << std::endl;
        const Scalar value = algo.integrate(f, canonicalSimplex)[0];
        fullprint << "n=" << n << " k=" << k << " value=" << std::setprecision(16) << value << std::endl;
        assert_almost_equal(value, ref[n]);
      }
    }
    // The integration over a uniform mesh of a symmetric integrand produces
    // many simplices with equal errors, so the selection of the simplices to
    // refine must break the ties with the simplex index: the result has to be
    // identical on every call and accurate.
    {
      IntervalMesher mesher(Indices(3, 8));
      Mesh mesh(mesher.build(Interval(Point(3, 0.0), Point(3, 1.0))));
      SymbolicFunction f(Description::BuildDefault(3, "x"), Description({"exp(-10.0*(x0+x1+x2))"}));
      const Scalar refValue = std::pow((1.0 - std::exp(-10.0)) / 10.0, 3);
      SimplicialCubature algo;
      algo.setRule(3);
      algo.setMaximumRelativeError(1.0e-7);
      const Point reference(algo.integrate(f, mesh));
      const UnsignedInteger callsAfterFirst = f.getCallsNumber();
      for (UnsignedInteger i = 0; i < 3; ++ i)
      {
        const Point value(algo.integrate(f, mesh));
        assert_almost_equal(value, reference);
        assert_almost_equal(value[0], refValue, 1.0e-6, 1.0e-8);
      }
      // The refinement must take place for the tie-breaking to be exercised
      const UnsignedInteger callsPerIntegration = f.getCallsNumber() - callsAfterFirst;
      if (!(callsPerIntegration > mesh.getSimplicesNumber() * 16))
        throw TestFailed(OSS() << "No refinement occurred, calls=" << callsPerIntegration);
      fullprint << "tie-breaking calls=" << f.getCallsNumber() << " value=" << std::setprecision(16) << reference[0] << std::endl;
    }
    // The batched evaluation must give the same result for any block size
    {
      IntervalMesher mesher(Indices(3, 8));
      Mesh mesh(mesher.build(Interval(Point(3, 0.0), Point(3, 1.0))));
      SymbolicFunction f(Description::BuildDefault(3, "x"), Description({"exp(-10.0*(x0+x1+x2))"}));
      const Scalar refValue = std::pow((1.0 - std::exp(-10.0)) / 10.0, 3);
      SimplicialCubature algo;
      algo.setRule(3);
      algo.setMaximumRelativeError(1.0e-6);
      Point reference;
      const UnsignedInteger savedBlockSize = ResourceMap::GetAsUnsignedInteger("SimplicialCubature-EvaluationBlockSize");
      try
      {
        for (const UnsignedInteger blockSize : {2048u, 7u, 1u})
        {
          ResourceMap::SetAsUnsignedInteger("SimplicialCubature-EvaluationBlockSize", blockSize);
          const Point value(algo.integrate(f, mesh));
          if (reference.getDimension() == 0) reference = value;
          assert_almost_equal(value[0], refValue, 1.0e-6, 1.0e-8);
          assert_almost_equal(value, reference, 1.0e-12, 0.0);
        }
      }
      catch (...)
      {
        ResourceMap::SetAsUnsignedInteger("SimplicialCubature-EvaluationBlockSize", savedBlockSize);
        throw;
      }
      ResourceMap::SetAsUnsignedInteger("SimplicialCubature-EvaluationBlockSize", savedBlockSize);
    }
    // The invalid inputs must be rejected
    {
      SimplicialCubature algo;
      Sample vertices(1, 2);
      vertices.add(Point({1.0, 0.0}));
      vertices.add(Point({0.0, 1.0}));
      Indices indices(3);
      indices[0] = 0;
      indices[1] = 1;
      indices[2] = 2;
      Mesh mesh(vertices, IndicesCollection(Collection<Indices>(1, indices)));
      bool caught = false;
      try
      {
        algo.integrate(SymbolicFunction(Description(1, "x"), Description(1, "x")), mesh);
      }
      catch (const InvalidArgumentException &)
      {
        caught = true;
      }
      assert_almost_equal(caught ? 1.0 : 0.0, 1.0);
      // An empty mesh of the right dimension integrates to zero
      Mesh emptyMesh(Sample(0, 2), IndicesCollection(0, 0));
      const Point value(algo.integrate(SymbolicFunction(Description::BuildDefault(2, "x"), Description({"x0", "x1"})), emptyMesh));
      assert_almost_equal(value, Point(2, 0.0));
      // A null block size must be rejected
      {
        const UnsignedInteger savedBlockSize = ResourceMap::GetAsUnsignedInteger("SimplicialCubature-EvaluationBlockSize");
        caught = false;
        try
        {
          ResourceMap::SetAsUnsignedInteger("SimplicialCubature-EvaluationBlockSize", 0);
          algo.integrate(SymbolicFunction(Description::BuildDefault(2, "x"), Description(1, "x0")), mesh);
        }
        catch (const InvalidArgumentException &)
        {
          caught = true;
        }
        catch (...)
        {
          ResourceMap::SetAsUnsignedInteger("SimplicialCubature-EvaluationBlockSize", savedBlockSize);
          throw;
        }
        ResourceMap::SetAsUnsignedInteger("SimplicialCubature-EvaluationBlockSize", savedBlockSize);
      }
      assert_almost_equal(caught ? 1.0 : 0.0, 1.0);
    }
  }
  catch (TestFailed & ex)
  {
    std::cerr << ex << std::endl;
    return ExitCode::Error;
  }
  return ExitCode::Success;
}
