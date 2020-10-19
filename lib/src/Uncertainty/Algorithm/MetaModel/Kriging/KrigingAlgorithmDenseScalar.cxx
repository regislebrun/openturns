//                                               -*- C++ -*-
/**
 *  @brief The class building gaussian process regression
 *
 *  Copyright 2005-2020 Airbus-EDF-IMACS-ONERA-Phimeca
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

#include "openturns/KrigingAlgorithmDenseScalar.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/LinearFunction.hxx"
#include "openturns/SpecFunc.hxx"
#include "openturns/KrigingEvaluation.hxx"
#include "openturns/KrigingGradient.hxx"
#include "openturns/CenteredFiniteDifferenceHessian.hxx"
#include "openturns/GeneralLinearModelResult.hxx"
#include "openturns/ComposedFunction.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(KrigingAlgorithmDenseScalar)

static const Factory<KrigingAlgorithmDenseScalar> Factory_KrigingAlgorithmDenseScalar;


/* Default constructor */
KrigingAlgorithmDenseScalar::KrigingAlgorithmDenseScalar()
  : MetaModelAlgorithm()
  , inputSample_(0, 0)
  , outputSample_(0, 0)
  , covarianceModel_()
  , glmAlgo_()
  , gamma_(0)
  , rho_(0)
  , result_()
  , covarianceCholeskyFactor_()
{
  // Nothing to do
}


/* Constructor */
KrigingAlgorithmDenseScalar::KrigingAlgorithmDenseScalar(const Sample & inputSample,
    const Sample & outputSample,
    const CovarianceModel & covarianceModel,
    const Basis & basis)
  : MetaModelAlgorithm()
  , inputSample_(inputSample)
  , outputSample_(outputSample)
  , covarianceModel_()
  , glmAlgo_(inputSample, outputSample, covarianceModel, basis, true)
  , gamma_(0)
  , rho_(0)
  , result_()
  , covarianceCholeskyFactor_()
{
  // Nothing to do
}


/* Virtual constructor */
KrigingAlgorithmDenseScalar * KrigingAlgorithmDenseScalar::clone() const
{
  return new KrigingAlgorithmDenseScalar(*this);
}

void KrigingAlgorithmDenseScalar::computeGamma()
{
  // Get cholesky factor & rho from glm
  LOGINFO("Solve L^t.gamma = rho");
  // Arguments are keepIntact=true, matrix_lower=true & solving_transposed=true
  gamma_ = covarianceCholeskyFactor_.getImplementation()->solveLinearSystemTri(rho_, true, true, true);
}

/* Perform regression */
void KrigingAlgorithmDenseScalar::run()
{
  LOGINFO("Launch GeneralLinearModelAlgorithmDenseScalar for the optimization");
  glmAlgo_.run();
  LOGINFO("End of GeneralLinearModelAlgorithm run");

  // Covariance coefficients are computed once, ever if optimiser is fixed
  rho_ = glmAlgo_.getRho();

  /* Method that returns the covariance factor */
  const GeneralLinearModelResult glmResult(glmAlgo_.getResult());
  covarianceCholeskyFactor_ = glmResult.getCholeskyFactor();
  LOGINFO("Compute the interpolation part");
  computeGamma();
  LOGINFO("Store the estimates");
  LOGINFO("Build the output meta-model");
  Function metaModel;
  // We use directly the collection of points
  const BasisCollection basis(glmResult.getBasisCollection());
  const CovarianceModel conditionalCovarianceModel(glmResult.getCovarianceModel());
  const Collection<Point> trendCoefficients(glmResult.getTrendCoefficients());
  const UnsignedInteger outputDimension = outputSample_.getDimension();
  Sample covarianceCoefficients(inputSample_.getSize(), outputDimension);
  covarianceCoefficients.getImplementation()->setData(gamma_);
  // Meta model definition
  metaModel.setEvaluation(new KrigingEvaluation(basis, inputSample_, conditionalCovarianceModel, trendCoefficients, covarianceCoefficients));
  metaModel.setGradient(new KrigingGradient(basis, inputSample_, conditionalCovarianceModel, trendCoefficients, covarianceCoefficients));
  metaModel.setHessian(new CenteredFiniteDifferenceHessian(ResourceMap::GetAsScalar( "CenteredFiniteDifferenceGradient-DefaultEpsilon" ), metaModel.getEvaluation()));

  // compute residual, relative error
  const Point outputVariance(outputSample_.computeVariance());
  
  const Sample mY(metaModel(inputSample_));
  //const Sample mY(outputSample_.getSize(), outputSample_.getDimension());
  const Point squaredResiduals((outputSample_ - mY).computeRawMoment(2));

  const UnsignedInteger size = inputSample_.getSize();
  Point residuals(outputDimension);
  Point relativeErrors(outputDimension);
  for (UnsignedInteger outputIndex = 0; outputIndex < outputDimension; ++ outputIndex)
  {
    residuals[outputIndex] = sqrt(squaredResiduals[outputIndex] / size);
    relativeErrors[outputIndex] = squaredResiduals[outputIndex] / outputVariance[outputIndex];
  }
  result_ = KrigingResult(inputSample_, outputSample_, metaModel, residuals, relativeErrors, basis, trendCoefficients, conditionalCovarianceModel, covarianceCoefficients, covarianceCholeskyFactor_);
}


/* Update regression */
void KrigingAlgorithmDenseScalar::update(const Sample & inputSample,
    const Sample & outputSample)
{
  glmAlgo_.run();
  const Sample oldInputSample(glmAlgo_.getInputSample());
  // Discretize the covariance model on the new points
  const CovarianceModel conditionalCovarianceModel(glmResult.getCovarianceModel());
  const UnsignedInteger oldSize = oldInputSample.getSize();
  const UnsignedInteger size = inputSample.getSize();
  const Matrix M12(oldSize, size);
  for (UnsignedInteger j = 0; j < size; ++j)
    for (UnsignedInteger i = 0; i < oldSize; ++i)
      M12(i, j) = conditionalCovarianceModel.computeAsScalar(inputSample[i], inputSample[j]);
  const CovarianceMatrix M22(conditionalCovarianceModel.discretize(inputSample));
  Matrix L12Transpose(covarianceCholeskyFactor_.solve(M12));
  Matrix L22(CovarianceMatrix((M22 - L12Transpose.computeGram()).getImplementation()).computeCholeski());

  /* Method that returns the covariance factor */
  const GeneralLinearModelResult glmResult(glmAlgo_.getResult());
  covarianceCholeskyFactor_ = glmResult.getCholeskyFactor();
  LOGINFO("Compute the interpolation part");
  computeGamma();
  LOGINFO("Store the estimates");
  LOGINFO("Build the output meta-model");
  Function metaModel;
  // We use directly the collection of points
  const BasisCollection basis(glmResult.getBasisCollection());
  const CovarianceModel conditionalCovarianceModel(glmResult.getCovarianceModel());
  const Collection<Point> trendCoefficients(glmResult.getTrendCoefficients());
  const UnsignedInteger outputDimension = outputSample_.getDimension();
  Sample covarianceCoefficients(inputSample_.getSize(), outputDimension);
  covarianceCoefficients.getImplementation()->setData(gamma_);
  // Meta model definition
  metaModel.setEvaluation(new KrigingEvaluation(basis, inputSample_, conditionalCovarianceModel, trendCoefficients, covarianceCoefficients));
  metaModel.setGradient(new KrigingGradient(basis, inputSample_, conditionalCovarianceModel, trendCoefficients, covarianceCoefficients));
  metaModel.setHessian(new CenteredFiniteDifferenceHessian(ResourceMap::GetAsScalar( "CenteredFiniteDifferenceGradient-DefaultEpsilon" ), metaModel.getEvaluation()));

  // compute residual, relative error
  const Point outputVariance(outputSample_.computeVariance());
  const Sample mY(metaModel(inputSample_));
  //const Sample mY(outputSample_.getSize(), outputSample_.getDimension());
  const Point squaredResiduals((outputSample_ - mY).computeRawMoment(2));

  const UnsignedInteger size = inputSample_.getSize();
  Point residuals(outputDimension);
  Point relativeErrors(outputDimension);
  for (UnsignedInteger outputIndex = 0; outputIndex < outputDimension; ++ outputIndex)
  {
    residuals[outputIndex] = sqrt(squaredResiduals[outputIndex] / size);
    relativeErrors[outputIndex] = squaredResiduals[outputIndex] / outputVariance[outputIndex];
  }
  result_ = KrigingResult(inputSample_, outputSample_, metaModel, residuals, relativeErrors, basis, trendCoefficients, conditionalCovarianceModel, covarianceCoefficients, covarianceCholeskyFactor_);
}


/* String converter */
String KrigingAlgorithmDenseScalar::__repr__() const
{
  return OSS() << "class=" << getClassName();
}


Sample KrigingAlgorithmDenseScalar::getInputSample() const
{
  return inputSample_;
}


Sample KrigingAlgorithmDenseScalar::getOutputSample() const
{
  return outputSample_;
}


KrigingResult KrigingAlgorithmDenseScalar::getResult()
{
  return result_;
}

/* Optimization solver accessor */
OptimizationAlgorithm KrigingAlgorithmDenseScalar::getOptimizationAlgorithm() const
{
  return glmAlgo_.getOptimizationAlgorithm();
}

void KrigingAlgorithmDenseScalar::setOptimizationAlgorithm(const OptimizationAlgorithm & solver)
{
  glmAlgo_.setOptimizationAlgorithm(solver);
}


/* Accessor to optimization bounds */
void KrigingAlgorithmDenseScalar::setOptimizationBounds(const Interval & optimizationBounds)
{
  glmAlgo_.setOptimizationBounds(optimizationBounds);
}

Interval KrigingAlgorithmDenseScalar::getOptimizationBounds() const
{
  return glmAlgo_.getOptimizationBounds();
}

/* Log-Likelihood function accessor */
Function KrigingAlgorithmDenseScalar::getReducedLogLikelihoodFunction()
{
  return glmAlgo_.getObjectiveFunction();
}

/* Optimize parameters flag accessor */
Bool KrigingAlgorithmDenseScalar::getOptimizeParameters() const
{
  return glmAlgo_.getOptimizeParameters();
}

void KrigingAlgorithmDenseScalar::setOptimizeParameters(const Bool optimizeParameters)
{
  glmAlgo_.setOptimizeParameters(optimizeParameters);
}

/* Observation noise accessor */
void KrigingAlgorithmDenseScalar::setNoise(const Point & noise)
{
  glmAlgo_.setNoise(noise);
}

Point KrigingAlgorithmDenseScalar::getNoise() const
{
  return glmAlgo_.getNoise();
}

/* Method save() stores the object through the StorageManager */
void KrigingAlgorithmDenseScalar::save(Advocate & adv) const
{
  MetaModelAlgorithm::save(adv);
  adv.saveAttribute( "inputSample_", inputSample_ );
  adv.saveAttribute( "outputSample_", outputSample_ );
  adv.saveAttribute( "covarianceModel_", covarianceModel_ );
  adv.saveAttribute( "result_", result_ );
  adv.saveAttribute( "covarianceCholeskyFactor_", covarianceCholeskyFactor_ );
}


/* Method load() reloads the object from the StorageManager */
void KrigingAlgorithmDenseScalar::load(Advocate & adv)
{
  MetaModelAlgorithm::load(adv);
  adv.loadAttribute( "inputSample_", inputSample_ );
  adv.loadAttribute( "outputSample_", outputSample_ );
  adv.loadAttribute( "covarianceModel_", covarianceModel_ );
  adv.loadAttribute( "result_", result_ );
  adv.loadAttribute( "covarianceCholeskyFactor_", covarianceCholeskyFactor_ );
}

END_NAMESPACE_OPENTURNS
