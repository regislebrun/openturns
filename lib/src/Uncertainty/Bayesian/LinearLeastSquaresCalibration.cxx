//                                               -*- C++ -*-
/**
 *  @brief Default LinearLeastSquaresCalibration
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
#include "openturns/LinearLeastSquaresCalibration.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/Normal.hxx"
#include "openturns/LinearFunction.hxx"
#include "openturns/SpecFunc.hxx"

BEGIN_NAMESPACE_OPENTURNS



CLASSNAMEINIT(LinearLeastSquaresCalibration)

static const Factory<LinearLeastSquaresCalibration> Factory_LinearLeastSquaresCalibration;

/* Default constructor */
LinearLeastSquaresCalibration::LinearLeastSquaresCalibration()
  : CalibrationAlgorithmImplementation()
{
  // Nothing to do
}

/* Parameter constructor */
LinearLeastSquaresCalibration::LinearLeastSquaresCalibration(const Function & model,
    const Sample & inputObservations,
    const Sample & outputObservations,
    const Point & startingPoint,
    const String & methodName)
  : CalibrationAlgorithmImplementation(model, inputObservations, outputObservations, Normal(startingPoint, CovarianceMatrix((IdentityMatrix(startingPoint.getDimension()) * SpecFunc::MaxScalar).getImplementation())))
  , modelObservations_(0, 0)
  , gradientObservations_(0, 0)
  , methodName_(methodName)
{
  // Check the input
  const UnsignedInteger parameterDimension = startingPoint.getDimension();
  if (model.getParameterDimension() != parameterDimension) throw InvalidArgumentException(HERE) << "Error: expected a model of parameter dimension=" << parameterDimension << ", got parameter dimension=" << model.getParameterDimension();
  const UnsignedInteger inputDimension = inputObservations_.getDimension();
  if (model.getInputDimension() != inputDimension) throw InvalidArgumentException(HERE) << "Error: expected a model of input dimension=" << inputDimension << ", got input dimension=" << model.getInputDimension();
  const UnsignedInteger outputDimension = outputObservations.getDimension();
  if (model.getOutputDimension() != outputDimension) throw InvalidArgumentException(HERE) << "Error: expected a model of output dimension=" << outputDimension << ", got output dimension=" << model.getOutputDimension();
  const UnsignedInteger size = outputObservations.getSize();
  if (outputObservations.getSize() != size) throw InvalidArgumentException(HERE) << "Error: expected an output sample of size=" << size << ", got size=" << outputObservations.getSize();
}

/* Parameter constructor */
LinearLeastSquaresCalibration::LinearLeastSquaresCalibration(const Sample & modelObservations,
    const Matrix & gradientObservations,
    const Sample & outputObservations,
    const Point & startingPoint,
    const String & methodName)
  : CalibrationAlgorithmImplementation(Function(), Sample(), outputObservations, Normal(startingPoint, CovarianceMatrix((IdentityMatrix(startingPoint.getDimension()) * SpecFunc::MaxScalar).getImplementation())))
  , modelObservations_(modelObservations)
  , gradientObservations_(gradientObservations)
  , methodName_(methodName)
{
  // Check the input
  const UnsignedInteger parameterDimension = startingPoint.getDimension();
  const UnsignedInteger outputDimension = outputObservations.getDimension();
  const UnsignedInteger size = outputObservations.getSize();
  if (gradientObservations.getNbColumns() != parameterDimension) throw InvalidArgumentException(HERE) << "Error: expected a gradient parameter of columns number=" << parameterDimension << ", got columns number=" << gradientObservations.getNbColumns();
  if (gradientObservations.getNbRows() != size * outputDimension) throw InvalidArgumentException(HERE) << "Error: expected a gradient parameter of rows number=" << size * outputDimension << ", got rows number=" << gradientObservations.getNbRows();
  if (gradientObservations.getNbColumns() != parameterDimension) throw InvalidArgumentException(HERE) << "Error: expected an observations gradient of column dimension=" << parameterDimension << ", got column dimension=" << gradientObservations.getNbColumns();
}

/* Performs the actual computation. Must be overloaded by the actual calibration algorithm */
void LinearLeastSquaresCalibration::run()
{
  if (getModel().getEvaluation().getImplementation()->isActualImplementation())
  {
    // Compute the linearization
    Function parametrizedModel(getModel());
    parametrizedModel.setParameter(getParameterPrior().getMean());
    // Flatten everything related to the model evaluations over the input observations
    const UnsignedInteger parameterDimension = getParameterPrior().getDimension();
    const UnsignedInteger outputDimension = getOutputObservations().getDimension();
    const UnsignedInteger size = getOutputObservations().getSize();
    modelObservations_ = parametrizedModel(getInputObservations());
    gradientObservations_ = MatrixImplementation(parameterDimension, size * outputDimension);
    UnsignedInteger shift = 0;
    UnsignedInteger skip = parameterDimension * outputDimension;
    for (UnsignedInteger i = 0; i < size; ++i)
    {
      const Matrix parameterGradient(parametrizedModel.parameterGradient(getInputObservations()[i]));
      std::copy(parameterGradient.getImplementation()->begin(), parameterGradient.getImplementation()->end(), gradientObservations_.getImplementation()->begin() + shift);
      shift += skip;
    }
    gradientObservations_ = gradientObservations_.transpose();
    parameterPrior_.setDescription(getModel().getParameterDescription());
  }

  const Point deltaY(modelObservations_.getImplementation()->getData() - outputObservations_.getImplementation()->getData());
  LeastSquaresMethod method(LeastSquaresMethod::Build(methodName_, gradientObservations_));
  const Point deltaTheta(method.solve(deltaY));
  for (UnsignedInteger i = 0; i < deltaTheta.getDimension(); ++ i)
    if (!SpecFunc::IsNormal(deltaTheta[i])) throw InvalidArgumentException(HERE) << "The calibration problem is not identifiable";

  const Point thetaStar(getStartingPoint() - deltaTheta);
  const Point r(deltaY - gradientObservations_ * deltaTheta);
  Distribution error;
  Distribution parameterPosterior(ComputePosteriorAndErrorDistribution(thetaStar, r, method.getGramInverse(), outputObservations_.getDimension(), error));
  parameterPosterior.setDescription(parameterPrior_.getDescription());
  const LinearFunction residualFunction(getStartingPoint(), deltaY, gradientObservations_);
  result_ = CalibrationResult(parameterPrior_, parameterPosterior, thetaStar, error, inputObservations_, outputObservations_, residualFunction, false);
}

Distribution LinearLeastSquaresCalibration::ComputePosteriorAndErrorDistribution(const Point & thetaStar,
    const Point & r,
    const SquareMatrix & gramInverse,
    const UnsignedInteger outputDimension,
    Distribution & error)
{
  const Scalar varianceError = r.normSquare() / (r.getDimension() - thetaStar.getDimension());
  CovarianceMatrix covarianceThetaStar;
  const Scalar epsilon = ResourceMap::GetAsScalar("LinearLeastSquaresCalibration-Regularization");
  covarianceThetaStar = CovarianceMatrix((gramInverse * varianceError).getImplementation());
  if (epsilon > 0.0)
  {
    const Scalar shift = epsilon * covarianceThetaStar.computeLargestEigenValueModule();
    for (UnsignedInteger i = 0; i < covarianceThetaStar.getDimension(); ++i)
      covarianceThetaStar(i, i) += shift;
  }
  Normal parameterPosterior;
  try
  {
    parameterPosterior = Normal(thetaStar, covarianceThetaStar);
  }
  catch (...)
  {
    throw InternalException(HERE) << "Error: the covariance of the posterior distribution is not definite positive. The problem may be not identifiable. Try to increase the \"LinearLeastSquaresCalibration-Regularization\" key in ResourceMap";
  }
  if (outputDimension > 0)
  {
    try
    {
      error = Normal(Point(outputDimension), (IdentityMatrix(outputDimension) * varianceError).getImplementation());
    }
    catch (...)
    {
      error = Normal(Point(outputDimension), (IdentityMatrix(outputDimension) * SpecFunc::MaxScalar).getImplementation());
    }
  } // outputDimension > 0
  return parameterPosterior;
}

/* Model observations accessor */
Sample LinearLeastSquaresCalibration::getModelObservations() const
{
  return modelObservations_;
}


/* Model gradient wrt the parameter accessor */
Matrix LinearLeastSquaresCalibration::getGradientObservations() const
{
  return gradientObservations_;
}

/* StartingPoint accessor */
Point LinearLeastSquaresCalibration::getStartingPoint() const
{
  // The startingPoint is stored in the prior distribution, which is a flat Normal distribution
  return getParameterPrior().getMean();
}

/* Least squares method name accessor */
String LinearLeastSquaresCalibration::getMethodName() const
{
  return methodName_;
}

/* String converter */
String LinearLeastSquaresCalibration::__repr__() const
{
  return OSS() << "class=" << LinearLeastSquaresCalibration::GetClassName()
         << " name=" << getName();
}


LinearLeastSquaresCalibration * LinearLeastSquaresCalibration::clone() const
{
  return new LinearLeastSquaresCalibration(*this);
}


/* Method save() stores the object through the StorageManager */
void LinearLeastSquaresCalibration::save(Advocate & adv) const
{
  PersistentObject::save(adv);
  adv.saveAttribute("modelObservations_", modelObservations_);
  adv.saveAttribute("gradientObservations_", gradientObservations_);
  adv.saveAttribute("methodName_", methodName_);
}

/* Method load() reloads the object from the StorageManager */
void LinearLeastSquaresCalibration::load(Advocate & adv)
{
  PersistentObject::load(adv);
  adv.loadAttribute("modelObservations_", modelObservations_);
  adv.loadAttribute("gradientObservations_", gradientObservations_);
  adv.loadAttribute("methodName_", methodName_);
}


END_NAMESPACE_OPENTURNS
