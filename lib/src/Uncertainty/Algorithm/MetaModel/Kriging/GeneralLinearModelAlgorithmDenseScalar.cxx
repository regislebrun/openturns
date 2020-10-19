//                                               -*- C++ -*-
/**
 *  @brief The class builds generalized linear models
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

#include "openturns/GeneralLinearModelAlgorithmDenseScalar.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/Log.hxx"
#include "openturns/SpecFunc.hxx"
#include "openturns/LinearFunction.hxx"
#include "openturns/NonCenteredFiniteDifferenceGradient.hxx"
#include "openturns/TNC.hxx"
#ifdef OPENTURNS_HAVE_ANALYTICAL_PARSER
#include "openturns/SymbolicFunction.hxx"
#else
#include "openturns/DatabaseFunction.hxx"
#endif
#include "openturns/IdentityFunction.hxx"
#include "openturns/ComposedFunction.hxx"
#include "openturns/LinearCombinationFunction.hxx"
#include "openturns/MemoizeFunction.hxx"
#include "openturns/DesignProxy.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(GeneralLinearModelAlgorithmDenseScalar)

static const Factory<GeneralLinearModelAlgorithmDenseScalar> Factory_GeneralLinearModelAlgorithmDenseScalar;

/* Default constructor */
GeneralLinearModelAlgorithmDenseScalar::GeneralLinearModelAlgorithmDenseScalar()
  : MetaModelAlgorithm()
  , inputSample_(0, 1) // 1 is to be consistent with the default covariance model
  , outputSample_(0, 1) // same
  , covarianceModel_()
  , reducedCovarianceModel_()
  , solver_()
  , optimizationBounds_()
  , beta_(0)
  , rho_(0)
  , F_(0, 0)
  , result_()
  , basis_()
  , covarianceCholeskyFactor_()
  , keepCholeskyFactor_(false)
  , hasRun_(false)
  , optimizeParameters_(true)
  , analyticalAmplitude_(false)
  , lastReducedLogLikelihood_(SpecFunc::LowestScalar)
{
  // Set the default covariance to adapt the active parameters of the covariance model
  setCovarianceModel(CovarianceModel());
  initializeDefaultOptimizationAlgorithm();
}

/* Parameters constructor */
GeneralLinearModelAlgorithmDenseScalar::GeneralLinearModelAlgorithmDenseScalar(const Sample & inputSample,
    const Sample & outputSample,
    const CovarianceModel & covarianceModel,
    const Bool keepCholeskyFactor)
  : MetaModelAlgorithm()
  , inputSample_(0, 0)
  , outputSample_(0, 0)
  , covarianceModel_()
  , reducedCovarianceModel_()
  , solver_()
  , optimizationBounds_()
  , beta_(0)
  , rho_(0)
  , F_(0, 0)
  , result_()
  , basis_()
  , covarianceCholeskyFactor_()
  , keepCholeskyFactor_(keepCholeskyFactor)
  , hasRun_(false)
  , optimizeParameters_(ResourceMap::GetAsBool("GeneralLinearModelAlgorithmDenseScalar-OptimizeParameters"))
  , analyticalAmplitude_(false)
  , lastReducedLogLikelihood_(SpecFunc::LowestScalar)
{
  // Set data
  setData(inputSample, outputSample);

  // If no basis then we suppose output sample centered
  checkYCentered(outputSample);

  // Set covariance model
  setCovarianceModel(covarianceModel);

  initializeDefaultOptimizationAlgorithm();
}

GeneralLinearModelAlgorithmDenseScalar::GeneralLinearModelAlgorithmDenseScalar(const Sample & inputSample,
    const Sample & outputSample,
    const CovarianceModel & covarianceModel,
    const Basis & basis,
    const Bool keepCholeskyFactor)
  : MetaModelAlgorithm()
  , inputSample_()
  , outputSample_()
  , covarianceModel_()
  , reducedCovarianceModel_()
  , solver_()
  , optimizationBounds_()
  , beta_(0)
  , rho_(0)
  , F_(0, 0)
  , result_()
  , basis_()
  , covarianceCholeskyFactor_()
  , keepCholeskyFactor_(keepCholeskyFactor)
  , method_(0)
  , hasRun_(false)
  , optimizeParameters_(ResourceMap::GetAsBool("GeneralLinearModelAlgorithmDenseScalar-OptimizeParameters"))
  , analyticalAmplitude_(false)
  , lastReducedLogLikelihood_(SpecFunc::LowestScalar)
{
  // Set data
  setData(inputSample, outputSample);

  // Set covariance model
  setCovarianceModel(covarianceModel);

  // Set basis
  setBasis(basis);

  initializeDefaultOptimizationAlgorithm();
}


/* set sample  method */
void GeneralLinearModelAlgorithmDenseScalar::setData(const Sample & inputSample,
    const Sample & outputSample)
{
  // Check the sample sizes
  if (inputSample.getSize() != outputSample.getSize())
    throw InvalidArgumentException(HERE) << "In GeneralLinearModelAlgorithmDenseScalar::GeneralLinearModelAlgorithmDenseScalar, input sample size=" << inputSample.getSize() << " does not match output sample size=" << outputSample.getSize();
  // Set samples
  inputSample_ = inputSample;
  outputSample_ = outputSample;
}

/* Covariance model accessors */
void GeneralLinearModelAlgorithmDenseScalar::setCovarianceModel(const CovarianceModel & covarianceModel)
{
  // Here we can store any modified version of the given covariance model wrt its parameters as it is mainly a parametric template

  // Here we focus only in the output dimension==1 case
  if (covarianceModel.getOutputDimension() != 1) throw InvalidArgumentException(HERE) << "Error: expected a covariance model with output dimension=1, got output dimension=" << covarianceModel.getOutputDimension();
  const UnsignedInteger inputDimension = inputSample_.getDimension();
  if (covarianceModel.getInputDimension() != inputDimension) throw InvalidArgumentException(HERE) << "Error: expected a covariance model with input dimension=1, got input dimension=" << covarianceModel.getInputDimension();

  // All the computation will be done on the reduced covariance model. We keep the initial covariance model (ie the one we just built) in order to reinitialize the reduced covariance model if some flags are changed after the creation of the algorithm.
  reducedCovarianceModel_ = covarianceModel_;
  // Now, adapt the model parameters.
  // First, check if the parameters have to be optimized. If not, remove all the active parameters.
  analyticalAmplitude_ = false;
  if (!optimizeParameters_) reducedCovarianceModel_.setActiveParameter(Indices());
  // Second, check if the amplitude parameter is unique and active
  else if (ResourceMap::GetAsBool("GeneralLinearModelAlgorithmDenseScalar-UseAnalyticalAmplitudeEstimate"))
  {
    const Description activeParametersDescription(reducedCovarianceModel_.getParameterDescription());
    // One of the active parameters must be called amplitude_0
    for (UnsignedInteger i = 0; i < activeParametersDescription.getSize(); ++i)
      if (activeParametersDescription[i] == "amplitude_0")
      {
        analyticalAmplitude_ = true;
        Indices newActiveParameters(reducedCovarianceModel_.getActiveParameter());
        newActiveParameters.erase(newActiveParameters.begin() + i);
        reducedCovarianceModel_.setActiveParameter(newActiveParameters);
        // Here we have to change the current value of the amplitude as it has
        // to be equal to 1 during the potential optimization step in order for
        // the analytical formula to be correct.
        // Now, the amplitude has disapear form the active parameters so it must
        // be updated using the amplitude accessor.
        reducedCovarianceModel_.setAmplitude(Point(1, 1.0));
        break;
      } // amplitude_0
  } // optimizeParameters_
  LOGINFO(OSS() << "final active parameters=" << reducedCovarianceModel_.getActiveParameter());
  // Define the bounds of the optimization problem
  const UnsignedInteger optimizationDimension = reducedCovarianceModel_.getParameter().getSize();
  if (optimizationDimension > 0)
  {
    const Scalar scaleFactor(ResourceMap::GetAsScalar( "GeneralLinearModelAlgorithmDenseScalar-DefaultOptimizationScaleFactor"));
    if (!(scaleFactor > 0))
      throw InvalidArgumentException(HERE) << "Scale factor set in ResourceMap is invalid. It should be a positive value. Here, scale = " << scaleFactor;
    const Point lowerBound(optimizationDimension, ResourceMap::GetAsScalar( "GeneralLinearModelAlgorithmDenseScalar-DefaultOptimizationLowerBound"));
    Point upperBound(optimizationDimension, ResourceMap::GetAsScalar( "GeneralLinearModelAlgorithmDenseScalar-DefaultOptimizationUpperBound"));
    // We can set scale parameter if these parameters are enabled.
    // check if scale is active
    const Indices activeParameters(reducedCovarianceModel_.getActiveParameter());
    Bool isScaleActive(true);
    for (UnsignedInteger k = 0; k < reducedCovarianceModel_.getScale().getSize(); ++k)
    {
      if (!activeParameters.contains(k))
        isScaleActive = false;
    }
    if (isScaleActive)
    {
      const Point inputSampleRange(inputSample_.getMax() - inputSample_.getMin());
      for (UnsignedInteger k = 0; k < reducedCovarianceModel_.getScale().getSize(); ++k) upperBound[k] = inputSampleRange[k] * scaleFactor;
    }
    LOGWARN(OSS() <<  "Warning! For coherency we set scale upper bounds = " << upperBound.__str__());

    optimizationBounds_ = Interval(lowerBound, upperBound);
  }
  else optimizationBounds_ = Interval();
}

CovarianceModel GeneralLinearModelAlgorithmDenseScalar::getCovarianceModel() const
{
  return covarianceModel_;
}

CovarianceModel GeneralLinearModelAlgorithmDenseScalar::getReducedCovarianceModel() const
{
  return reducedCovarianceModel_;
}

/* Set basis method */
void GeneralLinearModelAlgorithmDenseScalar::setBasis(const Basis & basis)
{
  // If basis given, then its size should be the same as the output dimension (each marginal of multibasis is a basis that will be used for trend of the corresponding marginal.
  if (basis.getSize() != outputSample_.getDimension())
    throw InvalidArgumentException(HERE) << "In GeneralLinearModelAlgorithmDenseScalar::GeneralLinearModelAlgorithmDenseScalar, output sample dimension=" << outputSample_.getDimension()  << " does not match multi-basis dimension=" << basis.getSize();
  if (basis.getOutputDimension() != 1) throw InvalidArgumentException(HERE) << "Error: expected a basis of output dimension=1, got output dimension=" << basis.getOutputDimension();
  basis_ = basis;
}

void GeneralLinearModelAlgorithmDenseScalar::checkYCentered(const Sample & Y)
{
  const Scalar meanEpsilon = ResourceMap::GetAsScalar("GeneralLinearModelAlgorithmDenseScalar-MeanEpsilon");
  const Point meanY(Y.computeMean());
  for (UnsignedInteger k = 0; k < meanY.getDimension(); ++k)
  {
    if (std::abs(meanY[k]) > meanEpsilon)
      LOGWARN(OSS() << "In GeneralLinearModelAlgorithmDenseScalar::GeneralLinearModelAlgorithmDenseScalar, basis is empty and output sample is not centered, mean=" << meanY);
  }
}

void GeneralLinearModelAlgorithmDenseScalar::initializeDefaultOptimizationAlgorithm()
{
  String solverName(ResourceMap::GetAsString("GeneralLinearModelAlgorithmDenseScalar-DefaultOptimizationAlgorithm"));
  solver_ = OptimizationAlgorithm::Build(solverName);
  TNC* tnc = dynamic_cast<TNC *>(solver_.getImplementation().get());
  if (tnc)
    tnc->setIgnoreFailure(true);
}

/* Virtual constructor */
GeneralLinearModelAlgorithmDenseScalar * GeneralLinearModelAlgorithmDenseScalar::clone() const
{
  return new GeneralLinearModelAlgorithmDenseScalar(*this);
}

/* Compute the design matrix */
void GeneralLinearModelAlgorithmDenseScalar::computeF()
{
  // Nothing to do if the design matrix has already been computed
  if (F_.getNbRows() != 0) return;
  LOGINFO("Compute the design matrix");
  // No early exit based on the sample/basis size as F_ must be initialized with the correct dimensions
  const UnsignedInteger sampleSize = inputSample_.getSize();
  const UnsignedInteger basisSize = basis.getSize();
  if (basisSize == 0) Matrix(sampleSize, basisSize);
  Indices indices(basis.getSize());
  indices.fill();
  designProxy_ = DesignProxy(inputSample_, basis.getSubBasis(indices));
  F_ = designProxy_.computeDesign(indices);
}


/* Perform regression
1) Compute the design matrix
2) Call the parameters optimization
  a) Compute the log-likelihood with the initial parameters. It is mandatory
     even if no parameter has to be optimized as this computation has many side
     effects such as:
     * computing the trend coefficients beta
     * computing the discretized covariance matrix Cholesky factor
  b) If the amplitude can be computed analytically from the other parameters:
     * set its value to 1
     * remove it from the list of parameters
  c) If some parameters remain, perform the optimization
  d) Deduce the associated value of the amplitude by the analytical formula if possible
3) Build the result:
  a) Extract the different parts of the trend
  b) Update the covariance model if needed
 */

void GeneralLinearModelAlgorithmDenseScalar::run()
{
  // Do not run again if already computed
  if (hasRun_) return;
  computeF();
  const UnsignedInteger outputDimension = outputSample_.getDimension();
  // optimization of likelihood function if provided
  // Here we call the optimizeReducedLogLikelihood() method even if the covariance
  // model has no active parameter, because:
  // + it can be due to the fact that the amplitude is obtained through an
  //   analytical formula and this situation is taken into account in
  //   maximizeReducedLogLikelihood()
  // + even if there is actually no parameter to optimize,
  //   maximizeReducedLogLikelihood() is the entry point to
  //   computeReducedLogLikelyhood() which has side effects on covariance
  //   discretization and factorization, and it computes beta_
  Scalar optimalLogLikelihood = maximizeReducedLogLikelihood();

  LOGINFO("Build the output meta-model");
  // The meta model is of type LinearCombination function
  Collection<Function> allFunctionsCollection(basis_);
  Function metaModel;

  if (allFunctionsCollection_.getSize() > 0)
  {
    // Care ! collection should be non empty
    metaModel = LinearCombinationFunction(allFunctionsCollectio, beta_);
  }
  else
  {
    // If no basis ==> zero function
#ifdef OPENTURNS_HAVE_ANALYTICAL_PARSER
    metaModel = SymbolicFunction(Description::BuildDefault(covarianceModel_.getInputDimension(), "x"), Description(covarianceModel_.getOutputDimension(), "0.0"));
#else
    metaModel = DatabaseFunction(Sample(1, reducedCovarianceModel_.getInputDimension()), Sample(1, reducedCovarianceModel_.getOutputDimension()));
#endif
  }

  // compute residual, relative error
  const Point outputVariance(outputSample_.computeVariance());
  const Sample mY(allFunctionsCollection_.getSize() == 0? Sample(inputSample_.getSize, 1) : metaModel(inputSample_));
  const Point squaredResiduals((outputSample_ - mY).computeRawMoment(2));

  Point residuals(outputDimension);
  Point relativeErrors(outputDimension);

  const UnsignedInteger size = inputSample_.getSize();
  for ( UnsignedInteger outputIndex = 0; outputIndex < outputDimension; ++ outputIndex )
  {
    residuals[outputIndex] = sqrt(squaredResiduals[outputIndex] / size);
    relativeErrors[outputIndex] = squaredResiduals[outputIndex] / outputVariance[outputIndex];
  }

  // return optimized covmodel with the original active parameters (see analyticalAmplitude_)
  CovarianceModel reducedCovarianceModelCopy(reducedCovarianceModel_);
  reducedCovarianceModelCopy.setActiveParameter(covarianceModel_.getActiveParameter());

  result_ = GeneralLinearModelResult(inputSample_, outputSample_, metaModel, residuals, relativeErrors, basisCollection_, trendCoefficients, reducedCovarianceModelCopy, optimalLogLikelihood);

  // The scaling is done there because it has to be done as soon as some optimization has been done, either numerically or through an analytical formula
  if (keepCholeskyFactor_)
  {
    if (analyticalAmplitude_)
    {
      const Scalar sigma = reducedCovarianceModel_.getAmplitude()[0];
      covarianceCholeskyFactor_ = covarianceCholeskyFactor_ * sigma;
    }
    result_.setCholeskyFactor(covarianceCholeskyFactor_);
  }
  else
    result_ = GeneralLinearModelResult(inputSample_, outputSample_, metaModel, residuals, relativeErrors, basisCollection_, trendCoefficients, reducedCovarianceModelCopy, optimalLogLikelihood);
  hasRun_ = true;
}

// Maximize the log-likelihood of the Gaussian process model wrt the observations
// If the covariance model has no active parameter, no numerical optimization
// is done. There are two cases:
// + no parameter has to be optimized, in which case a single call to
//   computeReducedLogLikelihood() is made in order to compute beta_ and to
//   factor the covariance matrix
// + the amplitude is the only covariance parameter to be estimated and it is
//   done thanks to an analytical formula
// The method returns the optimal log-likelihood (which is equal to the optimal
// reduced log-likelihood), the corresponding parameters being directly stored
// into the covariance model
Scalar GeneralLinearModelAlgorithmDenseScalar::maximizeReducedLogLikelihood()
{
  // initial guess
  Point initialParameters(reducedCovarianceModel_.getParameter());
  Indices initialActiveParameters(reducedCovarianceModel_.getActiveParameter());
  // We use the functional form of the log-likelihood computation to benefit from the cache mechanism
  Function reducedLogLikelihoodFunction(getObjectiveFunction());
  const Bool noNumericalOptimization = initialParameters.getSize() == 0;
  // Early exit if the parameters are known
  if (noNumericalOptimization)
  {
    // We only need to compute the log-likelihood function at the initial parameters in order to get the Cholesky factor and the trend coefficients
    const Scalar initialReducedLogLikelihood = reducedLogLikelihoodFunction(initialParameters)[0];
    LOGINFO("No covariance parameter to optimize");
    LOGINFO(OSS() << "initial parameters=" << initialParameters << ", log-likelihood=" << initialReducedLogLikelihood);
    return initialReducedLogLikelihood;
  }
  // At this point we have an optimization problem to solve
  // Define the optimization problem
  OptimizationProblem problem(reducedLogLikelihoodFunction);
  problem.setMinimization(false);
  problem.setBounds(optimizationBounds_);
  solver_.setStartingPoint(initialParameters);
  solver_.setProblem(problem);
  LOGINFO(OSS(false) << "Solve problem=" << problem << " using solver=" << solver_);
  solver_.run();
  const OptimizationAlgorithm::Result result(solver_.getResult());
  const Scalar optimalLogLikelihood = result.getOptimalValue()[0];
  const Point optimalParameters = result.getOptimalPoint();
  const UnsignedInteger evaluationNumber = result.getEvaluationNumber();
  // Check if the optimal value corresponds to the last computed value, in order to
  // see if the by-products (Cholesky factor etc) are correct
  if (lastReducedLogLikelihood_ != optimalLogLikelihood)
  {
    LOGDEBUG(OSS(false) << "Need to evaluate the objective function one more time because the last computed reduced log-likelihood value=" << lastReducedLogLikelihood_ << " is different from the optimal one=" << optimalLogLikelihood);
    (void) computeReducedLogLikelihood(optimalParameters);
  }
  // Final call to reducedLogLikelihoodFunction() in order to update the amplitude
  // No additional cost since the cache mechanism is activated
  LOGINFO(OSS() << evaluationNumber << " evaluations, optimized parameters=" << optimalParameters << ", log-likelihood=" << optimalLogLikelihood);

  return optimalLogLikelihood;
}

Point GeneralLinearModelAlgorithmDenseScalar::computeReducedLogLikelihood(const Point & parameters) const
{
  // Check that the parameters have a size compatible with the covariance model
  if (parameters.getSize() != reducedCovarianceModel_.getParameter().getSize())
    throw InvalidArgumentException(HERE) << "In GeneralLinearModelAlgorithmDenseScalar::computeReducedLogLikelihood, could not compute likelihood,"
                                         << " covariance model requires an argument of size " << reducedCovarianceModel_.getParameter().getSize()
                                         << " but here we got " << parameters.getSize();
  LOGDEBUG(OSS(false) << "Compute reduced log-likelihood for parameters=" << parameters);
  const Scalar constant = - SpecFunc::LOGSQRT2PI * static_cast<Scalar>(inputSample_.getSize()) * static_cast<Scalar>(outputSample_.getDimension());
  Scalar logDeterminant = 0.0;
  // If the amplitude is deduced from the other parameters, work with
  // the correlation function
  if (analyticalAmplitude_) reducedCovarianceModel_.setAmplitude(Point(1, 1.0));
  reducedCovarianceModel_.setParameter(parameters);
  // First, compute the log-determinant of the Cholesky factor of the covariance
  // matrix. As a by-product, also compute rho.
  logDeterminant = computeLapackLogDeterminantCholesky();
  // Compute the amplitude using an analytical formula if needed
  // and update the reduced log-likelihood.
  if (analyticalAmplitude_)
  {
    LOGDEBUG("Analytical amplitude");
    // J(\sigma)=-\log(\sqrt{\sigma^{2N}\det{R}})-(Y-M)^tR^{-1}(Y-M)/(2\sigma^2)
    //          =-N\log(\sigma)-\log(\det{R})/2-(Y-M)^tR^{-1}(Y-M)/(2\sigma^2)
    // dJ/d\sigma=-N/\sigma+(Y-M)^tR^{-1}(Y-M)/\sigma^3=0
    // \sigma=\sqrt{(Y-M)^tR^{-1}(Y-M)/N}
    const UnsignedInteger size = inputSample_.getSize();
    const Scalar sigma = std::sqrt(rho_.normSquare() / (ResourceMap::GetAsBool("GeneralLinearModelAlgorithmDenseScalar-UnbiasedVariance") ? size - beta_.getSize() : size));
    LOGDEBUG(OSS(false) << "sigma=" << sigma);
    reducedCovarianceModel_.setAmplitude(Point(1, sigma));
    logDeterminant += 2.0 * size * std::log(sigma);
    rho_ /= sigma;
    LOGDEBUG(OSS(false) << "rho_=" << rho_);
  } // analyticalAmplitude

  LOGDEBUG(OSS(false) << "log-determinant=" << logDeterminant << ", rho=" << rho_);
  const Scalar epsilon = rho_.normSquare();
  LOGDEBUG(OSS(false) << "epsilon=||rho||^2=" << epsilon);
  if (epsilon <= 0) lastReducedLogLikelihood_ = SpecFunc::LowestScalar;
  // For the general multidimensional case, we have to compute the general log-likelihood (ie including marginal variances)
  else lastReducedLogLikelihood_ = constant - 0.5 * (logDeterminant + epsilon);
  LOGINFO(OSS(false) << "Reduced log-likelihood=" << lastReducedLogLikelihood_);
  return Point(1, lastReducedLogLikelihood_);
}


Scalar GeneralLinearModelAlgorithmDenseScalar::computeLapackLogDeterminantCholesky() const
{
  // Using the hypothesis that parameters = scale & model writes : C(s,t) = diag(sigma) * R(s,t) * diag(sigma) with R a correlation function
  LOGDEBUG(OSS(false) << "Compute the LAPACK log-determinant of the Cholesky factor for covariance=" << reducedCovarianceModel_);

  LOGDEBUG("Discretize the covariance model");
  CovarianceMatrix C(reducedCovarianceModel_.discretize(inputSample_));
  LOGDEBUG(OSS(false) << "C=\n" << C);
  LOGDEBUG("Compute the Cholesky factor of the covariance matrix");
  Bool continuationCondition = true;
  const Scalar startingScaling = ResourceMap::GetAsScalar("GeneralLinearModelAlgorithmDenseScalar-StartingScaling");
  const Scalar maximalScaling = ResourceMap::GetAsScalar("GeneralLinearModelAlgorithmDenseScalar-MaximalScaling");
  Scalar cumulatedScaling = 0.0;
  Scalar scaling = startingScaling;
  while (continuationCondition && (cumulatedScaling < maximalScaling))
  {
    try
    {
      covarianceCholeskyFactor_ = C.computeCholesky();
      continuationCondition = false;
    }
    // If it has not yet been computed, compute it and store it
    catch (InternalException &)
    {
      cumulatedScaling += scaling ;
      // Unroll the regularization to optimize the computation
      for (UnsignedInteger i = 0; i < C.getDimension(); ++i) C(i, i) += scaling;
      scaling *= 2.0;
    }
  }
  if (scaling >= maximalScaling)
    throw InvalidArgumentException(HERE) << "In GeneralLinearModelAlgorithmDenseScalar::computeLapackLogDeterminantCholesky, could not compute the Cholesky factor."
                                         << " Scaling up to "  << cumulatedScaling << " was not enough";
  if (cumulatedScaling > 0.0)
    LOGWARN(OSS() <<  "Warning! Scaling up to "  << cumulatedScaling << " was needed in order to get an admissible covariance. ");
  LOGDEBUG(OSS(false) << "L=\n" << covarianceCholeskyFactor_);

  // y corresponds to output data
  const Point y(outputSample_.getImplementation()->getData());
  LOGDEBUG(OSS(false) << "y=" << y);
  // rho = L^{-1}y
  LOGDEBUG("Solve L.rho = y");
  rho_ = covarianceCholeskyFactor_.solveLinearSystem(y);
  LOGDEBUG(OSS(false) << "rho_=L^{-1}y=" << rho_);
  // If trend to estimate
  if (basisCollection_.getSize() > 0)
  {
    // Phi = L^{-1}F
    LOGDEBUG("Solve L.Phi = F");
    LOGDEBUG(OSS(false) << "F_=\n" << F_);
    Matrix Phi(covarianceCholeskyFactor_.solveLinearSystem(F_));
    LOGDEBUG(OSS(false) << "Phi=\n" << Phi);
    LOGDEBUG("Solve min_beta||Phi.beta - rho||^2");
    beta_ = Phi.solveLinearSystem(rho_);
    LOGDEBUG(OSS(false) << "beta_=" << beta_);
    LOGDEBUG("Update rho");
    rho_ -= Phi * beta_;
    LOGDEBUG(OSS(false) << "rho_=L^{-1}y-L^{-1}F.beta=" << rho_);
  }
  LOGDEBUG("Compute log(|det(L)|)=log(sqrt(|det(C)|))");
  Scalar logDetL = 0.0;
  for (UnsignedInteger i = 0; i < covarianceCholeskyFactor_.getDimension(); ++i )
  {
    const Scalar lii = covarianceCholeskyFactor_(i, i);
    if (lii <= 0.0) return SpecFunc::LowestScalar;
    logDetL += log(lii);
  }
  LOGDEBUG(OSS(false) << "logDetL=" << logDetL);
  return 2.0 * logDetL;
}

/* Optimization solver accessor */
OptimizationAlgorithm GeneralLinearModelAlgorithmDenseScalar::getOptimizationAlgorithm() const
{
  return solver_;
}

void GeneralLinearModelAlgorithmDenseScalar::setOptimizationAlgorithm(const OptimizationAlgorithm & solver)
{
  solver_ = solver;
  hasRun_ = false;
}

/* Optimize parameters flag accessor */
Bool GeneralLinearModelAlgorithmDenseScalar::getOptimizeParameters() const
{
  return optimizeParameters_;
}

void GeneralLinearModelAlgorithmDenseScalar::setOptimizeParameters(const Bool optimizeParameters)
{
  if (optimizeParameters != optimizeParameters_)
  {
    optimizeParameters_ = optimizeParameters;
    // Here we have to call setCovarianceModel() as it computes reducedCovarianceModel from covarianceModel_ in a way influenced by optimizeParameters_ flag.
    setCovarianceModel(covarianceModel_);
  }
}

/* Accessor to optimization bounds */
void GeneralLinearModelAlgorithmDenseScalar::setOptimizationBounds(const Interval & optimizationBounds)
{
  if (!(optimizationBounds.getDimension() == optimizationBounds_.getDimension())) throw InvalidArgumentException(HERE) << "Error: expected bounds of dimension=" << optimizationBounds_.getDimension() << ", got dimension=" << optimizationBounds.getDimension();
  optimizationBounds_ = optimizationBounds;
  hasRun_ = false;
}

Interval GeneralLinearModelAlgorithmDenseScalar::getOptimizationBounds() const
{
  return optimizationBounds_;
}

Point GeneralLinearModelAlgorithmDenseScalar::getRho() const
{
  return rho_;
}

/* String converter */
String GeneralLinearModelAlgorithmDenseScalar::__repr__() const
{
  OSS oss;
  oss << "class=" << getClassName()
      << ", inputSample=" << inputSample_
      << ", outputSample=" << outputSample_
      << ", basis=" << basisCollection_
      << ", covarianceModel=" << covarianceModel_
      << ", reducedCovarianceModel=" << reducedCovarianceModel_
      << ", solver=" << solver_
      << ", optimizeParameters=" << optimizeParameters_;
  return oss;
}


Sample GeneralLinearModelAlgorithmDenseScalar::getInputSample() const
{
  return inputSample_;
}


Sample GeneralLinearModelAlgorithmDenseScalar::getOutputSample() const
{
  return outputSample_;
}


GeneralLinearModelResult GeneralLinearModelAlgorithmDenseScalar::getResult()
{
  if (!hasRun_) run();
  return result_;
}


Function GeneralLinearModelAlgorithmDenseScalar::getObjectiveFunction()
{
  computeF();
  MemoizeFunction logLikelihood(ReducedLogLikelihoodEvaluation(*this));
  // Here we change the finite difference gradient for a non centered one in order to reduce the computational cost
  const Scalar finiteDifferenceEpsilon = ResourceMap::GetAsScalar( "NonCenteredFiniteDifferenceGradient-DefaultEpsilon" );
  logLikelihood.setGradient(NonCenteredFiniteDifferenceGradient(finiteDifferenceEpsilon, logLikelihood.getEvaluation()).clone());
  logLikelihood.enableCache();
  return logLikelihood;
}

void GeneralLinearModelAlgorithmDenseScalar::reset()
{
  // Reset elements for new computation
  // No need to update F_ as computeF /setBasisCollection are private
  // Same remark for setCovarianceModel & setData
  covarianceCholeskyFactor_ = TriangularMatrix();
  hasRun_ = false;
  lastReducedLogLikelihood_ = SpecFunc::LowestScalar;
}

/* Method save() stores the object through the StorageManager */
void GeneralLinearModelAlgorithmDenseScalar::save(Advocate & adv) const
{
  MetaModelAlgorithm::save(adv);
  adv.saveAttribute( "inputSample_", inputSample_ );
  adv.saveAttribute( "outputSample_", outputSample_ );
  adv.saveAttribute( "covarianceModel_", covarianceModel_ );
  adv.saveAttribute( "reducedCovarianceModel_", reducedCovarianceModel_ );
  adv.saveAttribute( "solver_", solver_ );
  adv.saveAttribute( "optimizationBounds_", optimizationBounds_ );
  adv.saveAttribute( "basis_", basis_ );
  adv.saveAttribute( "result_", result_ );
  adv.saveAttribute( "keepCholeskyFactor_", keepCholeskyFactor_ );
  adv.saveAttribute( "covarianceCholeskyFactor_", covarianceCholeskyFactor_ );
  adv.saveAttribute( "optimizeParameters_", optimizeParameters_ );
}


/* Method load() reloads the object from the StorageManager */
void GeneralLinearModelAlgorithmDenseScalar::load(Advocate & adv)
{
  MetaModelAlgorithm::load(adv);
  adv.loadAttribute( "inputSample_", inputSample_ );
  adv.loadAttribute( "outputSample_", outputSample_ );
  adv.loadAttribute( "covarianceModel_", covarianceModel_ );
  adv.loadAttribute( "reducedCovarianceModel_", reducedCovarianceModel_ );
  adv.loadAttribute( "solver_", solver_ );
  adv.loadAttribute( "optimizationBounds_", optimizationBounds_ );
  adv.loadAttribute( "basis_", basis_ );
  adv.loadAttribute( "result_", result_ );
  adv.loadAttribute( "method", method_ );
  adv.loadAttribute( "keepCholeskyFactor_", keepCholeskyFactor_ );
  adv.loadAttribute( "covarianceCholeskyFactor_", covarianceCholeskyFactor_ );
  adv.loadAttribute( "optimizeParameters_", optimizeParameters_ );
}

END_NAMESPACE_OPENTURNS
