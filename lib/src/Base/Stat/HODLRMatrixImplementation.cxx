//                                               -*- C++ -*-
/**
 *  @file  HODLRMatrixImplementation.cxx
 *  @brief HODLR compressed matrix implementation
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
 */
#include "openturns/HODLRMatrixImplementation.hxx"
#include "openturns/Log.hxx"
#include "openturns/OSS.hxx"

BEGIN_NAMESPACE_OPENTURNS

/** Evaluator that scales another evaluator by a constant factor. */
class HODLRScalingEvaluator
  : public HODLREntryEvaluator
{
public:
  HODLRScalingEvaluator(Pointer<const HODLREntryEvaluator> inner, Scalar alpha)
    : inner_(inner)
    , alpha_(alpha)
  {
  }

  HODLRScalingEvaluator* clone() const override
  {
    return new HODLRScalingEvaluator(inner_, alpha_);
  }

  Scalar operator()(UnsignedInteger i, UnsignedInteger j) const override
  {
    return alpha_ * (*inner_)(i, j);
  }

  UnsignedInteger getSize() const override
  {
    return inner_->getSize();
  }

private:
  Pointer<const HODLREntryEvaluator> inner_;
  Scalar alpha_;
};

CLASSNAMEINIT(HODLRMatrixImplementation)

HODLRMatrixImplementation::HODLRMatrixImplementation()
  : PersistentObject()
  , n_(0)
  , p_node_()
  , isFactorized_(false)
  , isCholesky_(false)
  , symmetric_(true)
  , logDet_(0.0)
  , shiftAccumulated_(0.0)
{
}

HODLRMatrixImplementation::HODLRMatrixImplementation(const HODLRMatrixImplementation& other)
  : PersistentObject(other)
  , n_(other.n_)
  , p_node_()
  , isFactorized_(false)
  , isCholesky_(false)
  , symmetric_(other.symmetric_)
  , logDet_(0.0)
  , parameters_(other.parameters_)
  , diagonal_(other.diagonal_)
  , shiftAccumulated_(0.0)
  , p_evaluator_(other.p_evaluator_)
{
}

HODLRMatrixImplementation* HODLRMatrixImplementation::clone() const
{
  return new HODLRMatrixImplementation(*this);
}

HODLRMatrixImplementation::~HODLRMatrixImplementation()
{
  // Nothing to do
}

HODLRMatrixImplementation& HODLRMatrixImplementation::operator=(const HODLRMatrixImplementation& other)
{
  if (this != &other)
  {
    PersistentObject::operator=(other);
    p_node_.reset();
    n_ = other.n_;
    isFactorized_ = false;
    isCholesky_ = false;
    symmetric_ = other.symmetric_;
    logDet_ = 0.0;
    parameters_ = other.parameters_;
    diagonal_ = other.diagonal_;
    shiftAccumulated_ = 0.0;
    p_evaluator_ = other.p_evaluator_;
    rng_ = std::mt19937();
  }
  return *this;
}

UnsignedInteger HODLRMatrixImplementation::getNbRows() const
{
  return n_;
}

UnsignedInteger HODLRMatrixImplementation::getNbColumns() const
{
  return n_;
}

const HODLRMatrixParameters& HODLRMatrixImplementation::getParameters() const
{
  return parameters_;
}

void HODLRMatrixImplementation::assemble(const HODLRRealAssemblyFunction& f, char symmetry)
{
  assemble(f, parameters_, symmetry);
}

void HODLRMatrixImplementation::assemble(const HODLRRealAssemblyFunction& f,
                                         const HODLRMatrixParameters& parameters,
                                         char symmetry)
{
  if (n_ == 0)
    throw InvalidArgumentException(HERE) << "Empty HODLRMatrix";

  // Validate the symmetry flag against the value given at build() time
  switch (symmetry)
  {
    case 'L':
    case 'l':
    case 'U':
    case 'u':
      if (!symmetric_)
        throw InvalidArgumentException(HERE) << "HODLRMatrix was built with symmetric=false but assembled with symmetry '" << symmetry << "'";
      break;
    case 'N':
    case 'n':
    case 'F':
    case 'f':
      if (symmetric_)
        throw InvalidArgumentException(HERE) << "HODLRMatrix was built with symmetric=true but assembled with symmetry '" << symmetry << "'";
      break;
    default:
      throw InvalidArgumentException(HERE) << "Invalid symmetry flag '" << symmetry << "', must be one of 'L', 'U', 'N' or 'F'";
  }

  parameters_ = parameters;

  p_evaluator_.reset(f.clone());

  // Build diagonal: the assembly function evaluates all entries including diagonal
  diagonal_ = Point(n_);
  for (UnsignedInteger i = 0; i < n_; ++i)
    diagonal_[i] = (*p_evaluator_)(i, i);

  // Build the HODLR tree
  p_node_.reset();
  isFactorized_ = false;

  rng_.seed(42);
  p_node_ = new HODLRNode(p_evaluator_, &diagonal_[0], 0, n_,
                        parameters.getMinLeafSize(),
                        parameters.getMaxRank(),
                        parameters.getAssemblyEpsilon(),
                        rng_);

  LOGDEBUG(OSS() << "HODLRMatrixImplementation::assemble done, n=" << n_);
}

void HODLRMatrixImplementation::rebuild()
{
  p_node_.reset();
  isFactorized_ = false;
  p_node_ = new HODLRNode(p_evaluator_, &diagonal_[0], 0, n_,
                          parameters_.getMinLeafSize(),
                          parameters_.getMaxRank(),
                          parameters_.getAssemblyEpsilon(),
                          rng_);
}

void HODLRMatrixImplementation::factorize(const String& method)
{
  if (n_ == 0)
    throw InvalidArgumentException(HERE) << "Empty HODLRMatrix";

  if (!p_node_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix not assembled";

  // Regularization loop: try factorization with increasing shift on failure.
  // Follows the same pattern as HMatrixImplementation::factorize().
  const Scalar regularizationEpsilon = ResourceMap::GetAsScalar("HODLRMatrix-RegularizationEpsilon");
  const UnsignedInteger maxIterations = ResourceMap::GetAsUnsignedInteger("HODLRMatrix-FactorizationIterations");
  if (maxIterations == 0)
    throw InvalidArgumentException(HERE) << "HODLRMatrix-FactorizationIterations must be positive";

  Scalar regularizationShift = shiftAccumulated_;
  for (UnsignedInteger iteration = 0; iteration < maxIterations; ++iteration)
  {
    p_node_->setShift(regularizationShift);
    try
    {
      if (method == "LLT" || method == "LLt")
      {
        isCholesky_ = true;
        p_node_->computeCholesky();
      }
      else
      {
        isCholesky_ = false;
        p_node_->compute();
      }
      shiftAccumulated_ = regularizationShift;
      break;
    }
    catch (const InternalException&)
    {
      LOGDEBUG(OSS() << "HODLRMatrix::factorize attempt " << iteration
                << " failed, regularizationShift=" << regularizationShift);
      if (regularizationShift == 0.0)
        regularizationShift = regularizationEpsilon;
      else
        regularizationShift *= 2.0;
      if (iteration == maxIterations - 1)
        throw InternalException(HERE) << "HODLRMatrix::factorize failed after "
                                      << maxIterations << " regularization attempts, last shift=" << regularizationShift;
      rebuild();
    }
  }

  logDet_ = p_node_->getLogDeterminant();
  isFactorized_ = true;
  LOGDEBUG(OSS() << "HODLRMatrix::factorize(" << method << ") done, log_det=" << logDet_);
}

void HODLRMatrixImplementation::scale(Scalar alpha)
{
  if (n_ == 0)
    throw InvalidArgumentException(HERE) << "Empty HODLRMatrix";

  if (!p_evaluator_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix not assembled";

  for (UnsignedInteger i = 0; i < n_; ++i)
    diagonal_[i] *= alpha;

  p_evaluator_ = new HODLRScalingEvaluator(p_evaluator_, alpha);
  rebuild();
}

void HODLRMatrixImplementation::gemv(char trans, Scalar alpha, const Point& x, Scalar beta, Point& y) const
{
  if (x.getDimension() != n_)
    throw InvalidArgumentException(HERE) << "x dimension mismatch";
  if (y.getDimension() != n_)
    throw InvalidArgumentException(HERE) << "y dimension mismatch";
  if (trans != 'N' && trans != 'n')
    throw InvalidArgumentException(HERE) << "trans must be 'N' (no transpose)";

  if (!p_node_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix not assembled";
  if (!isFactorized_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix not factorized";

  Matrix xmat(n_, 1);
  Matrix ymat(n_, 1);
  for (UnsignedInteger i = 0; i < n_; ++i)
  {
    xmat(i, 0) = x[i];
    ymat(i, 0) = 0.0;
  }

  if (isCholesky_)
    p_node_->applyFactor(ymat, xmat);
  else
    p_node_->apply(ymat, xmat);

  for (UnsignedInteger i = 0; i < n_; ++i)
    y[i] = alpha * ymat(i, 0) + beta * y[i];
}

void HODLRMatrixImplementation::addIdentity(Scalar alpha)
{
  shiftAccumulated_ += alpha;
  for (UnsignedInteger i = 0; i < n_; ++i)
    diagonal_[i] += alpha;
  if (p_node_)
    p_node_->setShift(shiftAccumulated_);
  isFactorized_ = false;
}

void HODLRMatrixImplementation::applyNugget()
{
  const Scalar nugget = ResourceMap::GetAsScalar("HODLRMatrix-Nugget");
  if (nugget <= 0.0)
    return;
  Scalar meanDiagonal = 0.0;
  for (UnsignedInteger i = 0; i < n_; ++i)
    meanDiagonal += diagonal_[i];
  meanDiagonal /= n_;
  addIdentity(nugget * meanDiagonal);
}

Scalar HODLRMatrixImplementation::norm() const
{
  throw NotYetImplementedException(HERE) << "HODLRMatrixImplementation::norm";
}

Point HODLRMatrixImplementation::getDiagonal() const
{
  return diagonal_;
}

Point HODLRMatrixImplementation::solve(const Point& b, Bool /* trans */) const
{
  if (!p_node_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix not assembled";

  if (!isFactorized_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix not factorized";

  Matrix x(b.getDimension(), 1);
  for (UnsignedInteger i = 0; i < b.getDimension(); ++i)
    x(i, 0) = b[i];

  p_node_->solve(x);

  Point result(b.getDimension());
  for (UnsignedInteger i = 0; i < b.getDimension(); ++i)
    result[i] = x(i, 0);
  return result;
}

Matrix HODLRMatrixImplementation::solve(const Matrix& m, Bool /* trans */) const
{
  if (!p_node_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix not assembled";

  if (!isFactorized_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix not factorized";

  Matrix x(m);
  p_node_->solve(x);
  return x;
}

Point HODLRMatrixImplementation::solveLower(const Point& b, Bool trans) const
{
  if (!p_node_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix not assembled";

  if (!isFactorized_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix not factorized";

  if (!isCholesky_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix must be Cholesky-factorized to use solveLower";

  Matrix x(b.getDimension(), 1);
  for (UnsignedInteger i = 0; i < b.getDimension(); ++i)
    x(i, 0) = b[i];

  p_node_->solveLower(x, trans);

  Point result(b.getDimension());
  for (UnsignedInteger i = 0; i < b.getDimension(); ++i)
    result[i] = x(i, 0);
  return result;
}

Matrix HODLRMatrixImplementation::solveLower(const Matrix& m, Bool trans) const
{
  if (!p_node_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix not assembled";

  if (!isFactorized_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix not factorized";

  if (!isCholesky_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix must be Cholesky-factorized to use solveLower";

  Matrix x(m);
  p_node_->solveLower(x, trans);
  return x;
}

Scalar HODLRMatrixImplementation::logDeterminant() const
{
  if (!p_node_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix not assembled";
  if (!isFactorized_)
    throw InvalidArgumentException(HERE) << "HODLRMatrix not factorized";
  // Cholesky: logDet_ holds log|L| = 0.5*log|A| (sum of log diag over all leaves)
  // LU: logDet_ already holds log|A| (sum of log abs diag over all blocks)
  return isCholesky_ ? 2.0 * logDet_ : logDet_;
}

std::pair<size_t, size_t> HODLRMatrixImplementation::compressionRatio() const
{
  if (!p_node_)
    return std::make_pair(static_cast<size_t>(0), static_cast<size_t>(n_) * n_);
  const size_t compressed = p_node_->getNnz();
  const size_t uncompressed = static_cast<size_t>(n_) * n_;
  return std::make_pair(compressed, uncompressed);
}

std::pair<size_t, size_t> HODLRMatrixImplementation::fullrkRatio() const
{
  const size_t full = static_cast<size_t>(n_) * n_;
  const size_t lr = p_node_ ? p_node_->getTotalRank() : 0;
  return std::make_pair(full, lr);
}

void HODLRMatrixImplementation::dump(const String& name) const
{
  LOGDEBUG(OSS() << "HODLRMatrix::dump(" << name << ")");
}

String HODLRMatrixImplementation::__repr__() const
{
  OSS oss(true);
  oss << "class= " << HODLRMatrixImplementation::GetClassName()
      << ", n= " << n_
      << ", factorized= " << isFactorized_;
  return oss;
}

String HODLRMatrixImplementation::__str__(const String&) const
{
  OSS oss(false);
  oss << "class= " << HODLRMatrixImplementation::GetClassName();
  return oss;
}

// --- HODLRCovarianceAssemblyFunction ---

HODLRCovarianceAssemblyFunction::HODLRCovarianceAssemblyFunction(
    const CovarianceModel& covarianceModel,
    const Sample& vertices)
  : HODLRRealAssemblyFunction()
  , covarianceModel_(covarianceModel)
  , vertices_(vertices)
  , verticesBegin_(vertices.getImplementation()->data_begin())
  , inputDimension_(vertices.getDimension())
  , covarianceDimension_(covarianceModel.getOutputDimension())
  , size_(vertices.getSize())
{
}

Scalar HODLRCovarianceAssemblyFunction::operator()(UnsignedInteger i, UnsignedInteger j) const
{
  if (covarianceDimension_ == 1)
    return covarianceModel_.getImplementation()->computeAsScalar(
        verticesBegin_ + i * inputDimension_, verticesBegin_ + j * inputDimension_);

  const UnsignedInteger rowIndex = i / covarianceDimension_;
  const UnsignedInteger columnIndex = j / covarianceDimension_;
  const SquareMatrix localCovarianceMatrix(covarianceModel_(vertices_[rowIndex], vertices_[columnIndex]));
  const UnsignedInteger rowIndexLocal = i % covarianceDimension_;
  const UnsignedInteger columnIndexLocal = j % covarianceDimension_;
  return localCovarianceMatrix(rowIndexLocal, columnIndexLocal);
}

UnsignedInteger HODLRCovarianceAssemblyFunction::getSize() const
{
  return covarianceDimension_ * size_;
}

END_NAMESPACE_OPENTURNS
