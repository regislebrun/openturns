//                                               -*- C++ -*-
/**
 *  @brief The test file of class IsoProbabilisticTransformation for standard methods
 *
 *  Copyright 2005-2025 Airbus-EDF-IMACS-ONERA-Phimeca
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

using namespace OT;
using namespace OT::Test;

typedef Distribution::IsoProbabilisticTransformation IsoProbabilisticTransformation;
typedef Distribution::InverseIsoProbabilisticTransformation InverseIsoProbabilisticTransformation;

SymmetricTensor clean(SymmetricTensor in)
{
  UnsignedInteger rowDim = in.getNbRows();
  UnsignedInteger colDim = in.getNbColumns();
  UnsignedInteger sheetDim = in.getNbSheets();
  for(UnsignedInteger i = 0; i < rowDim; i++)
    for(UnsignedInteger j = 0; j < colDim; j++)
      for (UnsignedInteger k = 0; k < sheetDim; k++)
      {
        in(i, j, k) = 1.e-4 * round(1.e4 * in(i, j, k));
        if (std::abs(in(i, j, k)) < 1.e-4) in(i, j, k) = 0.0;
      }
  return in;
}

int main(int, char *[])
{
  TESTPREAMBLE;
  OStream fullprint(std::cout);

  try
  {
    // Create a collection of distribution
    JointDistribution::DistributionCollection aCollection;

    aCollection.add(Uniform(-1.0, 2.0));
    aCollection.add(Gamma(2.0, 2.0, 0.0));

    UnsignedInteger dim = aCollection.getSize();

    // Create a copula
    CorrelationMatrix RCopula(dim);
    for (UnsignedInteger i = 0; i < dim; i++)
    {
      for (UnsignedInteger j = 0; j < i; j++)
      {
        RCopula(i, j) = (i + j + 1.0) / (2.0 * dim);
      }
    }
    JointDistribution distribution(aCollection, NormalCopula(RCopula));
    // Test for sampling
    UnsignedInteger size = 10000;
    Sample sample = distribution.getSample( size );
    fullprint << "sample first=" << sample[0] << " last=" << sample[size - 1] << std::endl;
    // Should be close to [0.5, 1, 0.544439]
    fullprint << "sample mean=" << sample.computeMean() << std::endl;
    // Should be close to [0.75, 0, 0, 0, 0.5, 0, 0, 0, 0.739866] for the independent copula
    fullprint << "sample covariance=" << sample.computeCovariance() << std::endl;

    IsoProbabilisticTransformation transform(distribution.getIsoProbabilisticTransformation());
    fullprint << "isoprobabilistic transformation=" << transform << std::endl;
    Sample transformedSample(transform(sample));
    fullprint << "transformed sample first=" << transformedSample[0] << " last=" << transformedSample[size - 1] << std::endl;
    fullprint << "transformed sample mean=" << transformedSample.computeMean() << std::endl;
    fullprint << "transformed sample covariance=" << transformedSample.computeCovariance() << std::endl;

    // Test for evaluation
    InverseIsoProbabilisticTransformation inverseTransform(distribution.getInverseIsoProbabilisticTransformation());
    fullprint << "inverse isoprobabilistic transformation=" << inverseTransform << std::endl;
    Sample transformedBackSample(inverseTransform(transformedSample));
    fullprint << "transformed back sample first=" << transformedBackSample[0] << " last=" << transformedBackSample[size - 1] << std::endl;
    fullprint << "transformed back sample mean=" << transformedBackSample.computeMean() << std::endl;
    fullprint << "transformed back sample covariance=" << transformedBackSample.computeCovariance() << std::endl;
    Point point(dim, 0.75);
    fullprint << "point=" << point << std::endl;
    Point transformedPoint(transform(point));
    fullprint << "transform value at point        =" << transformedPoint << std::endl;
    fullprint << "transform gradient at point     =" << transform.gradient(point) << std::endl;
    fullprint << "transform gradient at point (FD)=" << CenteredFiniteDifferenceGradient(1.0e-5, transform.getEvaluation()).gradient(point) << std::endl;
    fullprint << "transform hessian at point      =" << std::setprecision(4) << clean(transform.hessian(point)) << std::endl;
    fullprint << "transform hessian at point (FD) =" << std::setprecision(4) << clean(CenteredFiniteDifferenceHessian(1.0e-4, transform.getEvaluation()).hessian(point)) << std::endl;
    fullprint << "inverse transform value at transformed point        =" << inverseTransform(transformedPoint) << std::endl;
    fullprint << "inverse transform gradient at transformed point     =" << inverseTransform.gradient(transformedPoint) << std::endl;
    fullprint << "inverse transform gradient at transformed point (FD)=" << CenteredFiniteDifferenceGradient(1.0e-5, inverseTransform.getEvaluation()).gradient(transformedPoint) << std::endl;
    fullprint << "inverse transform hessian at transformed point      =" << std::setprecision(4) << clean(inverseTransform.hessian(transformedPoint)) << std::endl;
    fullprint << "inverse transform hessian at transformed point (FD) =" << std::setprecision(4) << clean(CenteredFiniteDifferenceHessian(1.0e-4, inverseTransform.getEvaluation()).hessian(transformedPoint)) << std::endl;

    // Test for parameters
    fullprint << "parameters gradient at point=" << transform.parameterGradient(point) << std::endl;
    {
      // Validation using finite difference
      Scalar eps = 1.0e-5;
      Scalar factor = 1.0 / (2.0 * eps);
      Matrix gradient(5, 2);
      Point dTdp;
      {
        // dT/dp0
        JointDistribution::DistributionCollection coll(aCollection);
        coll[0] = Uniform(-1.0 + eps, 2.0);
        IsoProbabilisticTransformation left(JointDistribution(coll, NormalCopula(RCopula)).getIsoProbabilisticTransformation());
        coll[0] = Uniform(-1.0 - eps, 2.0);
        IsoProbabilisticTransformation right(JointDistribution(coll, NormalCopula(RCopula)).getIsoProbabilisticTransformation());
        dTdp = (left(point) - right(point)) * factor;
        gradient(0, 0) = dTdp[0];
        gradient(0, 1) = dTdp[1];
      }
      {
        // dT/dp1
        JointDistribution::DistributionCollection coll(aCollection);
        coll[0] = Uniform(-1.0, 2.0 + eps);
        IsoProbabilisticTransformation left(JointDistribution(coll, NormalCopula(RCopula)).getIsoProbabilisticTransformation());
        coll[0] = Uniform(-1.0, 2.0 - eps);
        IsoProbabilisticTransformation right(JointDistribution(coll, NormalCopula(RCopula)).getIsoProbabilisticTransformation());
        dTdp = (left(point) - right(point)) * factor;
        gradient(1, 0) = dTdp[0];
        gradient(1, 1) = dTdp[1];
      }
      {
        // dT/dp2
        JointDistribution::DistributionCollection coll(aCollection);
        coll[1] = Gamma(2.0 + eps, 2.0, 0.0);
        IsoProbabilisticTransformation left(JointDistribution(coll, NormalCopula(RCopula)).getIsoProbabilisticTransformation());
        coll[1] = Gamma(2.0 - eps, 2.0, 0.0);
        IsoProbabilisticTransformation right(JointDistribution(coll, NormalCopula(RCopula)).getIsoProbabilisticTransformation());
        dTdp = (left(point) - right(point)) * factor;
        gradient(2, 0) = dTdp[0];
        gradient(2, 1) = dTdp[1];
      }
      {
        // dT/dp3
        JointDistribution::DistributionCollection coll(aCollection);
        coll[1] = Gamma(2.0, 2.0 + eps, 0.0);
        IsoProbabilisticTransformation left(JointDistribution(coll, NormalCopula(RCopula)).getIsoProbabilisticTransformation());
        coll[1] = Gamma(2.0, 2.0 - eps, 0.0);
        IsoProbabilisticTransformation right(JointDistribution(coll, NormalCopula(RCopula)).getIsoProbabilisticTransformation());
        dTdp = (left(point) - right(point)) * factor;
        gradient(3, 0) = dTdp[0];
        gradient(3, 1) = dTdp[1];
      }
      {
        // dT/dp4
        JointDistribution::DistributionCollection coll(aCollection);
        coll[1] = Gamma(2.0, 2.0, 0.0 + eps);
        IsoProbabilisticTransformation left(JointDistribution(coll, NormalCopula(RCopula)).getIsoProbabilisticTransformation());
        coll[1] = Gamma(2.0, 2.0, 0.0 - eps);
        IsoProbabilisticTransformation right(JointDistribution(coll, NormalCopula(RCopula)).getIsoProbabilisticTransformation());
        dTdp = (left(point) - right(point)) * factor;
        gradient(4, 0) = dTdp[0];
        gradient(4, 1) = dTdp[1];
      }
      fullprint << "parameters gradient (FD)    =" << gradient << std::endl;
    }
    // Test for parameters
    fullprint << "(inverse) parameters gradient at point=" << inverseTransform.parameterGradient(point) << std::endl;
    {
      // Validation using finite difference
      Scalar eps = 1.0e-5;
      Scalar factor = 1.0 / (2.0 * eps);
      Matrix gradient(5, 2);
      Point dTdp;
      {
        // dT/dp0
        JointDistribution::DistributionCollection coll(aCollection);
        coll[0] = Uniform(-1.0 + eps, 2.0);
        IsoProbabilisticTransformation left(JointDistribution(coll, NormalCopula(RCopula)).getInverseIsoProbabilisticTransformation());
        coll[0] = Uniform(-1.0 - eps, 2.0);
        IsoProbabilisticTransformation right(JointDistribution(coll, NormalCopula(RCopula)).getInverseIsoProbabilisticTransformation());
        dTdp = (left(point) - right(point)) * factor;
        gradient(0, 0) = dTdp[0];
        gradient(0, 1) = dTdp[1];
      }
      {
        // dT/dp1
        JointDistribution::DistributionCollection coll(aCollection);
        coll[0] = Uniform(-1.0, 2.0 + eps);
        IsoProbabilisticTransformation left(JointDistribution(coll, NormalCopula(RCopula)).getInverseIsoProbabilisticTransformation());
        coll[0] = Uniform(-1.0, 2.0 - eps);
        IsoProbabilisticTransformation right(JointDistribution(coll, NormalCopula(RCopula)).getInverseIsoProbabilisticTransformation());
        dTdp = (left(point) - right(point)) * factor;
        gradient(1, 0) = dTdp[0];
        gradient(1, 1) = dTdp[1];
      }
      {
        // dT/dp2
        JointDistribution::DistributionCollection coll(aCollection);
        coll[1] = Gamma(2.0 + eps, 2.0, 0.0);
        IsoProbabilisticTransformation left(JointDistribution(coll, NormalCopula(RCopula)).getInverseIsoProbabilisticTransformation());
        coll[1] = Gamma(2.0 - eps, 2.0, 0.0);
        IsoProbabilisticTransformation right(JointDistribution(coll, NormalCopula(RCopula)).getInverseIsoProbabilisticTransformation());
        dTdp = (left(point) - right(point)) * factor;
        gradient(2, 0) = dTdp[0];
        gradient(2, 1) = dTdp[1];
      }
      {
        // dT/dp3
        JointDistribution::DistributionCollection coll(aCollection);
        coll[1] = Gamma(2.0, 2.0 + eps, 0.0);
        IsoProbabilisticTransformation left(JointDistribution(coll, NormalCopula(RCopula)).getInverseIsoProbabilisticTransformation());
        coll[1] = Gamma(2.0, 2.0 - eps, 0.0);
        IsoProbabilisticTransformation right(JointDistribution(coll, NormalCopula(RCopula)).getInverseIsoProbabilisticTransformation());
        dTdp = (left(point) - right(point)) * factor;
        gradient(3, 0) = dTdp[0];
        gradient(3, 1) = dTdp[1];
      }
      {
        // dT/dp4
        JointDistribution::DistributionCollection coll(aCollection);
        coll[1] = Gamma(2.0, 2.0, 0.0 + eps);
        IsoProbabilisticTransformation left(JointDistribution(coll, NormalCopula(RCopula)).getInverseIsoProbabilisticTransformation());
        coll[1] = Gamma(2.0, 2.0, 0.0 - eps);
        IsoProbabilisticTransformation right(JointDistribution(coll, NormalCopula(RCopula)).getInverseIsoProbabilisticTransformation());
        dTdp = (left(point) - right(point)) * factor;
        gradient(4, 0) = dTdp[0];
        gradient(4, 1) = dTdp[1];
      }
      fullprint << "(inverse) parameters gradient (FD)    =" << gradient << std::endl;
    }
  }
  catch (TestFailed & ex)
  {
    std::cerr << ex << std::endl;
    return ExitCode::Error;
  }


  return ExitCode::Success;
}
