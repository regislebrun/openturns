//                                               -*- C++ -*-
/**
 *  @file  HODLRCore.cxx
 *  @brief Implementation of the recursive HODLR tree node
 *
 *  Based on the algorithm from george (Dan Foreman-Mackey, MIT license).
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
#include "openturns/HODLRCore.hxx"
#include "openturns/Lapack.hxx"
#include "openturns/Log.hxx"
#include "openturns/OSS.hxx"
#include "openturns/ResourceMap.hxx"

#include <algorithm>
#include <numeric>
#include <cmath>

BEGIN_NAMESPACE_OPENTURNS

BEGIN_C_DECLS
void dgetrs_(const char* trans, int* n, int* nrhs, double* a, int* lda,
             int* ipiv, double* b, int* ldb, int* info, int* ltrans);
void dpotrs_(const char* uplo, int* n, int* nrhs, double* a, int* lda,
             double* b, int* ldb, int* info, int* luplo);

#ifdef OPENTURNS_HAVE_OPENBLAS
int goto_get_num_procs(void);
void openblas_set_num_threads(int num_threads);
#endif
END_C_DECLS

namespace {

void HODLRDgemm(const char* transa, const char* transb, int* m, int* n, int* k,
                double* alpha, double* a, int* lda, double* b, int* ldb,
                double* beta, double* c, int* ldc)
{
  int ltransa = 1;
  int ltransb = 1;
  dgemm_(const_cast<char*>(transa), const_cast<char*>(transb), m, n, k, alpha,
         a, lda, b, ldb, beta, c, ldc, &ltransa, &ltransb);
}

void HODLRDtrsm(const char* side, const char* uplo, const char* transa, const char* diag,
                int* m, int* n, double* alpha, double* a, int* lda,
                double* b, int* ldb)
{
  int lside = 1;
  int luplo = 1;
  int ltransa = 1;
  int ldiag = 1;
  dtrsm_(const_cast<char*>(side), const_cast<char*>(uplo), const_cast<char*>(transa),
         const_cast<char*>(diag), m, n, alpha, a, lda, b, ldb,
         &lside, &luplo, &ltransa, &ldiag);
}

void HODLRDtrmv(const char* uplo, const char* transa, const char* diag,
                int* n, double* a, int* lda, double* x, int* incx)
{
  int luplo = 1;
  int ltransa = 1;
  int ldiag = 1;
  dtrmv_(const_cast<char*>(uplo), const_cast<char*>(transa), const_cast<char*>(diag),
         n, a, lda, x, incx, &luplo, &ltransa, &ldiag);
}

void HODLRDpotrf(const char* uplo, int* n, double* a, int* lda, int* info)
{
  int luplo = 1;
  dpotrf_(const_cast<char*>(uplo), n, a, lda, info, &luplo);
}

void HODLRDgetrs(const char* transa, int* n, int* nrhs, double* a, int* lda,
                 int* ipiv, double* b, int* ldb, int* info)
{
  int ltransa = 1;
  dgetrs_(const_cast<char*>(transa), n, nrhs, a, lda, ipiv, b, ldb, info, &ltransa);
}

void HODLRDpotrs(const char* uplo, int* n, int* nrhs, double* a, int* lda,
                 double* b, int* ldb, int* info)
{
  int luplo = 1;
  dpotrs_(const_cast<char*>(uplo), n, nrhs, a, lda, b, ldb, info, &luplo);
}

}  // namespace

HODLRNode::HODLRNode(Pointer<const HODLREntryEvaluator> eval,
                     const Scalar* diag,
                     UnsignedInteger start,
                     UnsignedInteger size,
                     UnsignedInteger minLeafSize,
                     UnsignedInteger maxRank,
                     Scalar tolerance,
                     std::mt19937& rng,
                     SignedInteger direction,
                     HODLRNode* parent)
  : p_diag_(diag)
  , p_eval_(eval)
  , start_(start)
  , size_(size)
  , direction_(direction)
  , rank_(0)
  , maxRank_(maxRank)
  , minLeafSize_(minLeafSize)
  , denseThreshold_(0)
  , tolerance_(tolerance)
  , p_rng_(&rng)
  , isLeaf_(false)
  , isCholesky_(false)
  , logDet_(0.0)
  , shift_(0.0)
  , totalRank_(0)
  , numLeaves_(0)
  , p_parent_(parent)
  , p_child0_()
  , p_child1_()
  , U_(2)
  , V_(2)
  , Uorig_(2)
{
  denseThreshold_ = ResourceMap::GetAsUnsignedInteger("HODLRMatrix-DenseThreshold");

  const UnsignedInteger half = size_ / 2;

  if (half >= minLeafSize)
  {
    isLeaf_ = false;

    // Low-rank approximation of off-diagonal block A01:
    // A[child1, child0] = U_[1] * V_[0]^T
    rank_ = lowRankApprox(
        start_ + half, size_ - half,
        start_, half,
        tolerance, rng, U_[1], V_[0]);

    // Symmetric: A10 = A01^T = V_[0] * U_[1]^T, so:
    U_[0] = V_[0];
    V_[1] = U_[1];

    totalRank_ = rank_;

    p_child0_ = new HODLRNode(p_eval_, p_diag_, start_, half, minLeafSize, maxRank_, tolerance, rng, 0, this);
    p_child1_ = new HODLRNode(p_eval_, p_diag_, start_ + half, size_ - half, minLeafSize, maxRank_, tolerance, rng, 1, this);
  }
  else
  {
    isLeaf_ = true;
  }
}

HODLRNode::~HODLRNode()
{
  // Nothing to do
}

void HODLRNode::setShift(Scalar shift)
{
  shift_ = shift;
  if (p_child0_) p_child0_->setShift(shift);
  if (p_child1_) p_child1_->setShift(shift);
}

HODLRBlasGuard::HODLRBlasGuard()
{
#ifdef OPENTURNS_HAVE_OPENBLAS
  savedNumThreads_ = goto_get_num_procs();
  openblas_set_num_threads(1);
#endif
  (void)savedNumThreads_;
}

HODLRBlasGuard::~HODLRBlasGuard()
{
#ifdef OPENTURNS_HAVE_OPENBLAS
  openblas_set_num_threads(savedNumThreads_);
#endif
}

UnsignedInteger HODLRNode::lowRankApprox(UnsignedInteger startRow, UnsignedInteger nRows,
                                         UnsignedInteger startCol, UnsignedInteger nCols,
                                         Scalar tol, std::mt19937& rng,
                                         Matrix& Uout, Matrix& Vout)
{
  const UnsignedInteger maxRankBound = std::min(nRows, nCols);
  const UnsignedInteger maxRank = (maxRank_ > 0) ? std::min(maxRankBound, maxRank_) : maxRankBound;

  std::vector<Scalar> Udata(nRows * maxRank, 0.0);
  std::vector<Scalar> Vdata(nCols * maxRank, 0.0);

  std::vector<UnsignedInteger> index(nRows);
  std::iota(index.begin(), index.end(), static_cast<UnsignedInteger>(0));
  std::shuffle(index.begin(), index.end(), rng);

  UnsignedInteger rank = 0;
  Scalar norm = 0.0;
  const Scalar tol2 = tol * tol;

  while (true)
  {
    SignedInteger i = -1;
    UnsignedInteger j = 0;

    while (true)
    {
      if (index.empty())
      {
        // All rows tested, none gave a valid pivot. Fall back to a
        // deterministic rank-maxRank approximation using the first columns.
        if ((maxRank_ > 0) && (maxRank_ < maxRankBound))
          LOGWARN(OSS() << "HODLRMatrix: rank-starved block (" << nRows << "x" << nCols
                 << ") reached the max rank " << maxRank_ << " before the assembly tolerance; "
                 << "accuracy may be degraded. Set HODLRMatrix-MaxRank to 0 for adaptive (tolerance-driven) rank");
        if (rank >= maxRank) break;
        Uout = Matrix(nRows, maxRank);
        Vout = Matrix(nCols, maxRank);
        MatrixImplementation& UoutImpl = *Uout.getImplementation();
        MatrixImplementation& VoutImpl = *Vout.getImplementation();
        if (nCols <= nRows)
        {
          const UnsignedInteger maxCols = std::min(nCols, maxRank);
          for (UnsignedInteger m = 0; m < maxCols; ++m)
            VoutImpl[m + m * nCols] = 1.0;
          for (UnsignedInteger n = 0; n < nRows; ++n)
            for (UnsignedInteger m = 0; m < maxCols; ++m)
              UoutImpl[n + m * nRows] = (*p_eval_)(startRow + n, startCol + m);
        }
        else
        {
          const UnsignedInteger maxRows = std::min(nRows, maxRank);
          for (UnsignedInteger m = 0; m < maxRows; ++m)
            UoutImpl[m + m * nRows] = 1.0;
          for (UnsignedInteger n = 0; n < maxRows; ++n)
            for (UnsignedInteger m = 0; m < nCols; ++m)
              VoutImpl[m + n * nCols] = (*p_eval_)(startRow + n, startCol + m);
        }
        return maxRank;
      }

      std::uniform_int_distribution<UnsignedInteger> dist(0, index.size() - 1);
      const UnsignedInteger k = dist(rng);
      i = static_cast<SignedInteger>(index[k]);
      index[k] = index.back();
      index.pop_back();

      for (UnsignedInteger n = 0; n < nCols; ++n)
        Vdata[n + rank * nCols] = (*p_eval_)(startRow + i, startCol + n);

      for (UnsignedInteger r = 0; r < rank; ++r)
      {
        const Scalar uir = Udata[i + r * nRows];
        for (UnsignedInteger n = 0; n < nCols; ++n)
          Vdata[n + rank * nCols] -= uir * Vdata[n + r * nCols];
      }

      Scalar maxVal = std::abs(Vdata[rank * nCols]);
      j = 0;
      for (UnsignedInteger n = 1; n < nCols; ++n)
      {
        const Scalar v = std::abs(Vdata[n + rank * nCols]);
        if (v > maxVal)
        {
          maxVal = v;
          j = n;
        }
      }

      if (maxVal >= 1e-14) break;
    }

    if (i < 0) break;

    const Scalar scale = 1.0 / Vdata[j + rank * nCols];
    for (UnsignedInteger n = 0; n < nCols; ++n)
      Vdata[n + rank * nCols] *= scale;

    for (UnsignedInteger n = 0; n < nRows; ++n)
      Udata[n + rank * nRows] = (*p_eval_)(startRow + n, startCol + j);

    for (UnsignedInteger r = 0; r < rank; ++r)
    {
      const Scalar vjr = Vdata[j + r * nCols];
      for (UnsignedInteger n = 0; n < nRows; ++n)
        Udata[n + rank * nRows] -= vjr * Udata[n + r * nRows];
    }

    ++rank;
    if (rank >= maxRank) break;

    Scalar uNorm2 = 0.0;
    Scalar vNorm2 = 0.0;
    for (UnsignedInteger n = 0; n < nRows; ++n)
      uNorm2 += Udata[n + (rank - 1) * nRows] * Udata[n + (rank - 1) * nRows];
    for (UnsignedInteger n = 0; n < nCols; ++n)
      vNorm2 += Vdata[n + (rank - 1) * nCols] * Vdata[n + (rank - 1) * nCols];
    Scalar rowcolNorm = uNorm2 * vNorm2;

    if (rowcolNorm < tol2 * norm) break;

    norm += rowcolNorm;
    if (rank > 1)
    {
      for (UnsignedInteger r = 0; r < rank - 1; ++r)
      {
        Scalar dotU = 0.0, dotV = 0.0;
        for (UnsignedInteger n = 0; n < nRows; ++n)
          dotU += Udata[n + r * nRows] * Udata[n + (rank - 1) * nRows];
        for (UnsignedInteger n = 0; n < nCols; ++n)
          dotV += Vdata[n + r * nCols] * Vdata[n + (rank - 1) * nCols];
        norm += 2.0 * dotU * dotV;
      }
    }
  }

  Uout = Matrix(nRows, rank);
  Vout = Matrix(nCols, rank);
  {
    MatrixImplementation& UoutImpl = *Uout.getImplementation();
    MatrixImplementation& VoutImpl = *Vout.getImplementation();
    std::copy(Udata.begin(), Udata.begin() + nRows * rank, &UoutImpl[0]);
    std::copy(Vdata.begin(), Vdata.begin() + nCols * rank, &VoutImpl[0]);
  }

  if ((maxRank_ > 0) && (rank >= maxRank_) && (rank < maxRankBound))
    LOGWARN(OSS() << "HODLRMatrix: rank-starved block (" << nRows << "x" << nCols
           << ") reached the max rank " << maxRank_ << " before the assembly tolerance; "
           << "accuracy may be degraded. Set HODLRMatrix-MaxRank to 0 for adaptive (tolerance-driven) rank");

  return rank;
}

void HODLRNode::compute()
{
  HODLRBlasGuard blasGuard;
  logDet_ = 0.0;

  if (!isLeaf_)
  {
    Uorig_[0] = U_[0];
    Uorig_[1] = U_[1];
    p_child0_->compute();
    p_child1_->compute();
    logDet_ = p_child0_->getLogDeterminant() + p_child1_->getLogDeterminant();
    totalRank_ = rank_ + p_child0_->getTotalRank() + p_child1_->getTotalRank();
    numLeaves_ = p_child0_->getNumLeaves() + p_child1_->getNumLeaves();
  }

  factorize();

  // Accumulate log-determinant from own factorization
  if (isLeaf_)
  {
    const MatrixImplementation& Sfact = *Sfactor_.getImplementation();
    for (UnsignedInteger n = 0; n < size_; ++n)
      logDet_ += std::log(std::abs(Sfact[n + n * size_]));
  }
  else
  {
    const UnsignedInteger dim = 2 * rank_;
    const MatrixImplementation& Sfact = *Sfactor_.getImplementation();
    for (UnsignedInteger n = 0; n < dim; ++n)
      logDet_ += std::log(std::abs(Sfact[n + n * dim]));
  }

  // Propagate apply_inverse up the parent chain (george's approach)
  HODLRNode* node = p_parent_;
  UnsignedInteger start = start_;
  SignedInteger ind = direction_;
  while (node)
  {
    applyInverse(node->U_[ind], start);
    start = node->start_;
    ind = node->direction_;
    node = node->p_parent_;
  }
}

// --- HODLRCorrectedEvaluator factory implementations ---

Pointer<HODLRCorrectedEvaluator> HODLRCorrectedEvaluator::create(
    Pointer<const HODLREntryEvaluator> original,
    UnsignedInteger offset,
    UnsignedInteger size,
    const Matrix& U1,
    const Matrix& K,
    Scalar lambda)
{
  const UnsignedInteger nRows = U1.getNbRows();
  const UnsignedInteger rank = K.getNbRows();
  Matrix UK(nRows, rank);
  if (rank > 0)
  {
    int m = static_cast<int>(nRows);
    int k = static_cast<int>(rank);
    int l = static_cast<int>(rank);
    double one = 1.0, zero = 0.0;
    HODLRDgemm("N", "N", &m, &k, &l, &one,
           const_cast<double*>(&U1(0, 0)), &m,
           const_cast<double*>(&K(0, 0)), &l,
           &zero, &UK(0, 0), &m);
  }

  Correction corr;
  corr.offset = offset;
  corr.size = size;
  corr.rank = rank;
  corr.U1 = Matrix(nRows, rank);
  if (rank > 0)
    std::copy(&U1(0, 0), &U1(0, 0) + nRows * rank, &corr.U1(0, 0));
  corr.UK = UK;
  corr.lambda = lambda;

  std::vector<Correction> corrections;
  corrections.push_back(corr);

  return new HODLRCorrectedEvaluator(original, offset, size, corrections);
}

Pointer<HODLRCorrectedEvaluator> HODLRCorrectedEvaluator::flatten(
    const HODLRCorrectedEvaluator& existing,
    const Matrix& newU1,
    const Matrix& newK,
    Scalar newLambda,
    UnsignedInteger newOffset,
    UnsignedInteger newSize)
{
  // Build UK = newU1 * newK via dgemm
  const UnsignedInteger nRows = newU1.getNbRows();
  const UnsignedInteger rank = newK.getNbRows();
  Matrix UK(nRows, rank);
  if (rank > 0)
  {
    int m = static_cast<int>(nRows);
    int k = static_cast<int>(rank);
    int l = static_cast<int>(rank);
    double one = 1.0, zero = 0.0;
    HODLRDgemm("N", "N", &m, &k, &l, &one,
           const_cast<double*>(&newU1(0, 0)), &m,
           const_cast<double*>(&newK(0, 0)), &l,
           &zero, &UK(0, 0), &m);
  }

  // Copy existing corrections
  std::vector<Correction> allCorrections = existing.corrections_;

  // Append new correction
  Correction newCorr;
  newCorr.offset = newOffset;
  newCorr.size = newSize;
  newCorr.rank = rank;
  newCorr.U1 = Matrix(nRows, rank);
  if (rank > 0)
    std::copy(&newU1(0, 0), &newU1(0, 0) + nRows * rank, &newCorr.U1(0, 0));
  newCorr.UK = UK;
  newCorr.lambda = newLambda;
  allCorrections.push_back(newCorr);

  return new HODLRCorrectedEvaluator(existing.getOriginal(),
                                     existing.getOffset(),
                                     existing.getSize(),
                                     allCorrections);
}

void HODLRNode::computeCholesky()
{
  HODLRBlasGuard blasGuard;
  isCholesky_ = true;

  if (isLeaf_)
  {
    factorizeLeafCholesky();
    logDet_ = 0.0;
    {
      const MatrixImplementation& Sfact = *Sfactor_.getImplementation();
      for (UnsignedInteger n = 0; n < size_; ++n)
        logDet_ += std::log(Sfact[n + n * size_]);
    }
    return;
  }

  Uorig_[0] = U_[0];
  Uorig_[1] = U_[1];

  // 1. Recursively compute L00 = chol(A00)
  p_child0_->computeCholesky();

  totalRank_ = rank_ + p_child0_->getTotalRank();
  numLeaves_ = p_child0_->getNumLeaves();

  const UnsignedInteger s0 = size_ / 2;
  const UnsignedInteger s1 = size_ - s0;

  if (rank_ > 0)
  {
    // 2. Compute W = L00^{-1} * V_[0]
    W_ = V_[0];
    p_child0_->applyInverseFactor(W_);

    // 3. Build K = W^T * W  (rank_ x rank_) via dgemm
    Matrix K(rank_, rank_);
    {
      int mK = static_cast<int>(rank_);
      int nK = static_cast<int>(rank_);
      int kK = static_cast<int>(s0);
      double one = 1.0, zero = 0.0;
      int ldW = static_cast<int>(s0);
      int ldK = static_cast<int>(rank_);
      HODLRDgemm("T", "N", &mK, &nK, &kK, &one, &W_(0, 0), &ldW, &W_(0, 0), &ldW, &zero, &K(0, 0), &ldK);
    }

    // 4. Schur complement: factorize A11' = A11 - U1 * K * U1^T
    // Start with no regularization; on first failure jump directly to a
    // meaningful value (1e-4) to avoid wasting attempts on tiny lambdas.
    // Then double geometrically, capped at maxLambda0.
    Scalar lambda = 0.0;
    const Bool useDenseFallback = (denseThreshold_ > 0) && (s1 <= denseThreshold_);
    const Scalar maxLambda0 = ResourceMap::GetAsScalar("HODLRMatrix-MaxRegularization");
    const Scalar lambdaFactor = ResourceMap::GetAsScalar("HODLRMatrix-RegularizationFactor");
    const UnsignedInteger maxAttempts = ResourceMap::GetAsUnsignedInteger("HODLRMatrix-RegularizationAttempts");
    if (maxAttempts == 0)
      throw InvalidArgumentException(HERE) << "HODLRMatrix-RegularizationAttempts must be positive";
    for (UnsignedInteger attempt = 0; attempt < maxAttempts; ++attempt)
    {
      if (useDenseFallback)
      {
        p_child1_ = new HODLRNode(p_eval_, p_diag_, start_ + s0, s1,
                                  s1 + 1, maxRank_, tolerance_, *p_rng_, 1, this);
        p_child1_->setShift(shift_);
        try
        {
          p_child1_->factorizeLeafCholeskyCorrected(K, U_[1], lambda);
          numLeaves_ += 1;
          break;
        }
        catch (const InternalException&)
        {
          if (lambda == 0.0)
            lambda = 1e-12;
          else if (lambda < maxLambda0)
            lambda = std::min(lambda * lambdaFactor, maxLambda0);
          else
            lambda *= lambdaFactor;
          if (attempt == maxAttempts - 1)
            throw InternalException(HERE) << "Dense Schur complement not SPD even with lambda=" << lambda;
        }
      }
      else
      {
        // Hierarchical Schur complement: rebuild child1 with FLATTENED corrected evaluator.
        Pointer<HODLRCorrectedEvaluator> correctedEval;
        const HODLRCorrectedEvaluator* existingEval =
            dynamic_cast<const HODLRCorrectedEvaluator*>(p_eval_.get());
        if (existingEval)
        {
          correctedEval = HODLRCorrectedEvaluator::flatten(
              *existingEval, U_[1], K, lambda, start_ + s0, s1);
        }
        else
        {
          correctedEval = HODLRCorrectedEvaluator::create(
              p_eval_, start_ + s0, s1, U_[1], K, lambda);
        }
        p_child1_ = new HODLRNode(correctedEval, p_diag_, start_ + s0, s1,
                                  minLeafSize_, maxRank_, tolerance_, *p_rng_, 1, this);
        p_child1_->setShift(shift_);
        try
        {
          p_child1_->computeCholesky();
          totalRank_ += p_child1_->getTotalRank();
          numLeaves_ += p_child1_->getNumLeaves();
          break;
        }
        catch (const InternalException&)
        {
          if (lambda == 0.0)
            lambda = 1e-12;
          else if (lambda < maxLambda0)
            lambda = std::min(lambda * lambdaFactor, maxLambda0);
          else
            lambda *= lambdaFactor;
          if (attempt == maxAttempts - 1)
            throw InternalException(HERE) << "Hierarchical Schur complement not SPD even with lambda=" << lambda;
        }
      }
    }
  }
  else
  {
    p_child1_->computeCholesky();
    totalRank_ += p_child1_->getTotalRank();
    numLeaves_ += p_child1_->getNumLeaves();
  }

  logDet_ = p_child0_->getLogDeterminant() + p_child1_->getLogDeterminant();

  // Build the S matrix for internal coupling (needed by parent applyInverse)
  if (!isLeaf_)
  {
    const UnsignedInteger dim = 2 * rank_;
    Sfactor_ = Matrix(dim, dim);
    MatrixImplementation& Sfact = *Sfactor_.getImplementation();
    for (UnsignedInteger i = 0; i < dim; ++i)
      for (UnsignedInteger j = 0; j < dim; ++j)
        Sfact[i + j * dim] = (i == j) ? 1.0 : 0.0;

    if (rank_ > 0)
    {
      const Scalar* V1_data = &(*V_[1].getImplementation())[0];
      const Scalar* U1_data = &(*U_[1].getImplementation())[0];
      const Scalar* V0_data = &(*V_[0].getImplementation())[0];
      const Scalar* U0_data = &(*U_[0].getImplementation())[0];
      {
        int m = static_cast<int>(rank_);
        int n = static_cast<int>(rank_);
        int k = static_cast<int>(s1);
        double one = 1.0;
        double zero = 0.0;
        int ld = static_cast<int>(dim);
        HODLRDgemm("T", "N", &m, &n, &k, &one,
               const_cast<double*>(V1_data), &k,
               const_cast<double*>(U1_data), &k,
               &zero, &Sfact[0 + rank_ * dim], &ld);
      }
      {
        int m = static_cast<int>(rank_);
        int n = static_cast<int>(rank_);
        int k = static_cast<int>(s0);
        double one = 1.0;
        double zero = 0.0;
        int ld = static_cast<int>(dim);
        HODLRDgemm("T", "N", &m, &n, &k, &one,
               const_cast<double*>(V0_data), &k,
               const_cast<double*>(U0_data), &k,
               &zero, &Sfact[rank_ + 0 * dim], &ld);
      }
    }

    ipiv_.resize(dim);
    int info = 0;
    int nn = static_cast<int>(dim);
    dgetrf_(&nn, &nn, &Sfact[0], &nn, ipiv_.data(), &info);
    if (info != 0)
      throw InternalException(HERE) << "LU factorization failed, info=" << info;
  }
}

void HODLRNode::factorize()
{
  if (isLeaf_)
  {
    // Build dense block: K(i,j) = kernel(start+i, start+j)
    Sfactor_ = Matrix(size_, size_);
    {
      MatrixImplementation& Sfact = *Sfactor_.getImplementation();
      for (UnsignedInteger i = 0; i < size_; ++i)
        for (UnsignedInteger j = 0; j < size_; ++j)
          Sfact[i + j * size_] = (*p_eval_)(start_ + i, start_ + j);
      if (shift_ != 0.0)
        for (UnsignedInteger i = 0; i < size_; ++i)
          Sfact[i + i * size_] += shift_;
    }
    leafMatrix_ = Matrix(size_, size_);
    {
      const Scalar* Sfact = &(*Sfactor_.getImplementation())[0];
      MatrixImplementation& leafData = *leafMatrix_.getImplementation();
      const UnsignedInteger leafSize = size_ * size_;
      std::copy(Sfact, Sfact + leafSize, &leafData[0]);
    }
    // LU factorize (in-place, leafMatrix_ preserves original for apply())
    ipiv_.resize(size_);
    int info = 0;
    int nn = static_cast<int>(size_);
    dgetrf_(&nn, &nn, &(*Sfactor_.getImplementation())[0], &nn, ipiv_.data(), &info);
    if (info != 0)
      throw InternalException(HERE) << "LU factorization failed, info=" << info;
  }
  else
  {
    // Build inner matrix S = I + [0, V1^T*U1; V0^T*U0, 0]
    const UnsignedInteger dim = 2 * rank_;

    Sfactor_ = Matrix(dim, dim);
    MatrixImplementation& Sfact = *Sfactor_.getImplementation();

    for (UnsignedInteger i = 0; i < dim; ++i)
      for (UnsignedInteger j = 0; j < dim; ++j)
        Sfact[i + j * dim] = (i == j) ? 1.0 : 0.0;

    if (rank_ > 0)
    {
      const UnsignedInteger s1 = size_ - size_ / 2;
      const Scalar* V1_data = &(*V_[1].getImplementation())[0];
      const Scalar* U1_data = &(*U_[1].getImplementation())[0];
      {
        int m = static_cast<int>(rank_);
        int n = static_cast<int>(rank_);
        int k = static_cast<int>(s1);
        double one = 1.0;
        double zero = 0.0;
        int ldS = static_cast<int>(dim);
        HODLRDgemm("T", "N", &m, &n, &k, &one,
               const_cast<double*>(V1_data), &k,
               const_cast<double*>(U1_data), &k,
               &zero, &Sfact[0 + rank_ * dim], &ldS);
      }

      const UnsignedInteger s0 = size_ / 2;
      const Scalar* V0_data = &(*V_[0].getImplementation())[0];
      const Scalar* U0_data = &(*U_[0].getImplementation())[0];
      {
        int m = static_cast<int>(rank_);
        int n = static_cast<int>(rank_);
        int k = static_cast<int>(s0);
        double one = 1.0;
        double zero = 0.0;
        int ldS = static_cast<int>(dim);
        HODLRDgemm("T", "N", &m, &n, &k, &one,
               const_cast<double*>(V0_data), &k,
               const_cast<double*>(U0_data), &k,
               &zero, &Sfact[rank_ + 0 * dim], &ldS);
      }
    }

    // LU factorize S
    ipiv_.resize(dim);
    int info = 0;
    int nn = static_cast<int>(dim);
    dgetrf_(&nn, &nn, &Sfact[0], &nn, ipiv_.data(), &info);
    if (info != 0)
      throw InternalException(HERE) << "LU factorization failed, info=" << info;
  }
}

void HODLRNode::factorizeLeafCholesky()
{
  if (!isLeaf_)
  {
    // Build S matrix (same as LU case -- non-symmetric, for Woodbury coupling)
    const UnsignedInteger dim = 2 * rank_;
    Sfactor_ = Matrix(dim, dim);
    MatrixImplementation& Sfact = *Sfactor_.getImplementation();
    for (UnsignedInteger i = 0; i < dim; ++i)
      for (UnsignedInteger j = 0; j < dim; ++j)
        Sfact[i + j * dim] = (i == j) ? 1.0 : 0.0;

    if (rank_ > 0)
    {
      const UnsignedInteger s1 = size_ - size_ / 2;
      const Scalar* V1_data = &(*V_[1].getImplementation())[0];
      const Scalar* U1_data = &(*U_[1].getImplementation())[0];
      {
        int m = static_cast<int>(rank_);
        int n = static_cast<int>(rank_);
        int k = static_cast<int>(s1);
        double one = 1.0;
        double zero = 0.0;
        int ldS = static_cast<int>(dim);
        HODLRDgemm("T", "N", &m, &n, &k, &one,
               const_cast<double*>(V1_data), &k,
               const_cast<double*>(U1_data), &k,
               &zero, &Sfact[0 + rank_ * dim], &ldS);
      }

      const UnsignedInteger s0 = size_ / 2;
      const Scalar* V0_data = &(*V_[0].getImplementation())[0];
      const Scalar* U0_data = &(*U_[0].getImplementation())[0];
      {
        int m = static_cast<int>(rank_);
        int n = static_cast<int>(rank_);
        int k = static_cast<int>(s0);
        double one = 1.0;
        double zero = 0.0;
        int ldS = static_cast<int>(dim);
        HODLRDgemm("T", "N", &m, &n, &k, &one,
               const_cast<double*>(V0_data), &k,
               const_cast<double*>(U0_data), &k,
               &zero, &Sfact[rank_ + 0 * dim], &ldS);
      }
    }

    // LU factorize S (S is not SPD, LU is correct for the Woodbury coupling)
    ipiv_.resize(dim);
    int info = 0;
    int nn = static_cast<int>(dim);
    dgetrf_(&nn, &nn, &Sfact[0], &nn, ipiv_.data(), &info);
    if (info != 0)
      throw InternalException(HERE) << "LU factorization failed, info=" << info;
    return;
  }

  // Leaf: assemble dense matrix
  Sfactor_ = Matrix(size_, size_);
  {
    MatrixImplementation& Sfact = *Sfactor_.getImplementation();

    // Check if evaluator has precomputed correction data (HODLRCorrectedEvaluator).
    // If so, use the ROOT original evaluator and apply all corrections via batch dgemm.
    const HODLRCorrectedEvaluator* corrEval =
        dynamic_cast<const HODLRCorrectedEvaluator*>(p_eval_.get());

    if (corrEval)
    {
      const Pointer<const HODLREntryEvaluator> original = corrEval->getOriginal();
      const auto& corrections = corrEval->getCorrections();
      const UnsignedInteger n = size_;
      const UnsignedInteger leafStart = start_;

      // Fill with root original evaluator entries
      for (UnsignedInteger i = 0; i < n; ++i)
        for (UnsignedInteger j = 0; j < n; ++j)
          Sfact[i + j * n] = (*original)(leafStart + i, leafStart + j);

      // Apply global shift to all diagonals (from addIdentity)
      if (shift_ != 0.0)
        for (UnsignedInteger i = 0; i < n; ++i)
          Sfact[i + i * n] += shift_;

      // Apply each correction via batch dgemm
      for (const auto& corr : corrections)
      {
        const SignedInteger localOffset =
            static_cast<SignedInteger>(corr.offset) - static_cast<SignedInteger>(leafStart);
        const SignedInteger localEnd =
            localOffset + static_cast<SignedInteger>(corr.size);
        if (localEnd <= 0) continue;   // correction entirely before leaf
        if (localOffset >= static_cast<SignedInteger>(n)) break;  // correction beyond leaf

        const UnsignedInteger i0 = static_cast<UnsignedInteger>(std::max(SignedInteger(0), localOffset));
        const UnsignedInteger i1 = static_cast<UnsignedInteger>(std::min(static_cast<SignedInteger>(n), localEnd));
        const UnsignedInteger localN = i1 - i0;

        if (localN == 0 || corr.rank == 0) continue;

        const UnsignedInteger ukRowOffset = i0 - static_cast<UnsignedInteger>(localOffset);

        // Sfact[i0:i1, i0:i1] -= UK[ukRowOffset:ukRowOffset+localN, :] * U1[ukRowOffset:ukRowOffset+localN, :]^T
        int m = static_cast<int>(localN);
        int k = static_cast<int>(corr.rank);
        int nMat = static_cast<int>(localN);
        double neg = -1.0, one = 1.0;
        int ldUK = static_cast<int>(corr.size);
        int ldU1 = static_cast<int>(corr.size);
        int ldS = static_cast<int>(n);
        HODLRDgemm("N", "T", &m, &nMat, &k, &neg,
               const_cast<double*>(&corr.UK(ukRowOffset, 0)), &ldUK,
               const_cast<double*>(&corr.U1(ukRowOffset, 0)), &ldU1,
               &one, &Sfact[i0 + i0 * n], &ldS);

        // Add lambda to diagonal within correction range
        if (corr.lambda != 0.0)
        {
          for (UnsignedInteger ii = i0; ii < i1; ++ii)
            Sfact[ii + ii * n] += corr.lambda;
        }
      }
    }
    else
    {
      for (UnsignedInteger i = 0; i < size_; ++i)
        for (UnsignedInteger j = 0; j < size_; ++j)
          Sfact[i + j * size_] = (*p_eval_)(start_ + i, start_ + j);
      if (shift_ != 0.0)
        for (UnsignedInteger i = 0; i < size_; ++i)
          Sfact[i + i * size_] += shift_;
    }
  }
  leafMatrix_ = Matrix(size_, size_);
  {
    const Scalar* Sfact = &(*Sfactor_.getImplementation())[0];
    MatrixImplementation& leafData = *leafMatrix_.getImplementation();
    const UnsignedInteger leafSize = size_ * size_;
    std::copy(Sfact, Sfact + leafSize, &leafData[0]);
  }

  int info = 0;
  int n = static_cast<int>(size_);
  {
    MatrixImplementation& Sfact2 = *Sfactor_.getImplementation();
    HODLRDpotrf("L", &n, &Sfact2[0], &n, &info);
  }
  if (info != 0)
    throw InternalException(HERE) << "Cholesky factorization failed, info=" << info;
}

void HODLRNode::factorizeLeafCholeskyCorrected(const Matrix& K, const Matrix& U1, Scalar lambda)
{
  // Directly assemble A11' = A11 - U1 * K * U1^T + lambda * I
  // Uses the ORIGINAL evaluator (p_eval_) for A11 entries, computes
  // the low-rank correction via dgemm for numerical stability.
  const UnsignedInteger n = size_;
  const UnsignedInteger rank = K.getNbRows();

  Sfactor_ = Matrix(n, n);
  MatrixImplementation& Sfact = *Sfactor_.getImplementation();

  // 1. Fill with original evaluator entries
  for (UnsignedInteger i = 0; i < n; ++i)
    for (UnsignedInteger j = 0; j < n; ++j)
      Sfact[i + j * n] = (*p_eval_)(start_ + i, start_ + j);

  // 2. Compute correction = U1 * K * U1^T via dgemm
  if (rank > 0)
  {
    // temp = U1 * K  (n x rank)
    Matrix temp(n, rank);
    {
      int m = static_cast<int>(n);
      int k = static_cast<int>(rank);
      int l = static_cast<int>(rank);
      double one = 1.0, zero = 0.0;
      HODLRDgemm("N", "N", &m, &k, &l, &one,
             const_cast<double*>(&U1(0, 0)), &m,
             const_cast<double*>(&K(0, 0)), &l,
             &zero, &temp(0, 0), &m);
    }
    // Sfact -= temp * U1^T  (n x n)
    {
      int m = static_cast<int>(n);
      int k = static_cast<int>(rank);
      int l = static_cast<int>(n);
      double neg = -1.0, one = 1.0;
      HODLRDgemm("N", "T", &m, &l, &k, &neg,
             &temp(0, 0), &m,
             const_cast<double*>(&U1(0, 0)), &m,
             &one, &Sfact[0], &m);
    }
  }

  // 3. Add lambda regularization to diagonal
  if (lambda != 0.0)
  {
    for (UnsignedInteger i = 0; i < n; ++i)
      Sfact[i + i * n] += lambda;
  }

  // 3b. Add global shift (from addIdentity) to diagonal
  if (shift_ != 0.0)
  {
    for (UnsignedInteger i = 0; i < n; ++i)
      Sfact[i + i * n] += shift_;
  }

  leafMatrix_ = Matrix(n, n);
  {
    MatrixImplementation& leafData = *leafMatrix_.getImplementation();
    const UnsignedInteger leafSize = n * n;
    std::copy(&Sfact[0], &Sfact[0] + leafSize, &leafData[0]);
  }

  // 4. Cholesky factorization
  int info = 0;
  int nn = static_cast<int>(n);
  HODLRDpotrf("L", &nn, &Sfact[0], &nn, &info);
  if (info != 0)
    throw InternalException(HERE) << "Cholesky factorization failed, info=" << info;

  // 5. Accumulate log-det contribution: sum(log(diag(L))) for this leaf
  logDet_ = 0.0;
  for (UnsignedInteger i = 0; i < n; ++i)
    logDet_ += std::log(Sfact[i + i * n]);
}

void HODLRNode::applyInverse(Matrix& x, UnsignedInteger start) const
{
  const UnsignedInteger nrhs = x.getNbColumns();
  const SignedInteger offset = static_cast<SignedInteger>(start_) - static_cast<SignedInteger>(start);

  if (isLeaf_)
  {
    if (isCholesky_)
    {
      applyInverseCholeskyLeaf(x, start);
    }
    else
    {
      applyInverseLeaf(x, start);
    }
    return;
  }

  const UnsignedInteger s0 = size_ / 2;
  const UnsignedInteger s1 = size_ - s0;

  if (rank_ == 0)
    return;

  const Scalar* V1_data = &(*V_[1].getImplementation())[0];
  const Scalar* V0_data = &(*V_[0].getImplementation())[0];
  const Scalar* U0_data = &(*U_[0].getImplementation())[0];
  const Scalar* U1_data = &(*U_[1].getImplementation())[0];

  // temp[0:rank_, :] = V1^T * x[s0:s0+s1, :]
  // temp[rank_:2*rank_, :] = V0^T * x[0:s0, :]
  Matrix temp(2 * rank_, nrhs, 0.0);
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger i = 0; i < rank_; ++i)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger r = 0; r < s1; ++r)
        sum += V1_data[r + i * s1] * x(offset + s0 + r, j);
      temp(i, j) = sum;
    }
    for (UnsignedInteger i = 0; i < rank_; ++i)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger r = 0; r < s0; ++r)
        sum += V0_data[r + i * s0] * x(offset + r, j);
      temp(rank_ + i, j) = sum;
    }
  }

  // Solve S * temp = temp
  if (rank_ > 0)
  {
    int n = static_cast<int>(2 * rank_);
    int nr = static_cast<int>(nrhs);
    int ldb = n;
    int info = 0;
    char trans = 'N';
    Matrix Scopy(Sfactor_);
    std::vector<int> ipivCopy(ipiv_.begin(), ipiv_.end());
    HODLRDgetrs(&trans, &n, &nr, &Scopy(0, 0), &n, ipivCopy.data(), &temp(0, 0), &ldb, &info);
    if (info != 0)
      throw InternalException(HERE) << "Solve failed, info=" << info;
  }

  // x[0:s0] -= U0 * temp[0:rank_]
  // x[s0:s0+s1] -= U1 * temp[rank_:2*rank_]
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger r = 0; r < s0; ++r)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger i = 0; i < rank_; ++i)
        sum += U0_data[r + i * s0] * temp(i, j);
      x(offset + r, j) -= sum;
    }
    for (UnsignedInteger r = 0; r < s1; ++r)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger i = 0; i < rank_; ++i)
        sum += U1_data[r + i * s1] * temp(rank_ + i, j);
      x(offset + s0 + r, j) -= sum;
    }
  }
}

void HODLRNode::applyInverseLeaf(Matrix& x, UnsignedInteger start) const
{
  const UnsignedInteger nrhs = x.getNbColumns();
  const SignedInteger offset = static_cast<SignedInteger>(start_) - static_cast<SignedInteger>(start);
  int n = static_cast<int>(size_);
  int nr = static_cast<int>(nrhs);
  int ldb = n;
  int info = 0;
  char trans = 'N';

  Matrix block(size_, nrhs);
  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < size_; ++i)
      block(i, j) = x(offset + i, j);

  Matrix Dcopy(Sfactor_);
  std::vector<int> ipivCopy(ipiv_.begin(), ipiv_.end());
  HODLRDgetrs(&trans, &n, &nr, &Dcopy(0, 0), &n, ipivCopy.data(), &block(0, 0), &ldb, &info);
  if (info != 0)
    throw InternalException(HERE) << "Solve failed, info=" << info;

  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < size_; ++i)
      x(offset + i, j) = block(i, j);
}

void HODLRNode::applyInverseCholeskyLeaf(Matrix& x, UnsignedInteger start) const
{
  const UnsignedInteger nrhs = x.getNbColumns();
  const SignedInteger offset = static_cast<SignedInteger>(start_) - static_cast<SignedInteger>(start);
  int n = static_cast<int>(size_);
  int nr = static_cast<int>(nrhs);
  int ldb = n;
  int info = 0;

  Matrix block(size_, nrhs);
  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < size_; ++i)
      block(i, j) = x(offset + i, j);

  // Solve L * L^T * block = block  =>  block = C^{-1} * block
  Matrix Lcopy(Sfactor_);
  HODLRDpotrs("L", &n, &nr, &Lcopy(0, 0), &n, &block(0, 0), &ldb, &info);
  if (info != 0)
    throw InternalException(HERE) << "Solve failed, info=" << info;

  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < size_; ++i)
      x(offset + i, j) = block(i, j);
}

void HODLRNode::applyInverseFactor(Matrix& x) const
{
  // Applies L^{-1} * x in-place where L is this node's HODLR Cholesky factor.
  // x must have size_ rows.
  const UnsignedInteger nrhs = x.getNbColumns();

  if (isLeaf_)
  {
    // Dense lower-triangular solve via dtrsm
    int n = static_cast<int>(size_);
    int nrhs_int = static_cast<int>(nrhs);
    double one = 1.0;
    HODLRDtrsm("L", "L", "N", "N", &n, &nrhs_int, &one,
           const_cast<double*>(&Sfactor_(0, 0)), &n, &x(0, 0), &n);
    return;
  }

  const UnsignedInteger s0 = size_ / 2;
  const UnsignedInteger s1 = size_ - s0;

  // Forward substitution with the 2x2 block structure:
  // [L00   0 ] * [y0] = [x0]
  // [L10  L11]   [y1]   [x1]
  //
  // 1. y0 = L00^{-1} * x0   (recursive on child0)
  // 2. y1 = L11^{-1} * (x1 - L10 * y0)
  //    where L10 = U_[1] * W_^T  (if isCholesky_ && rank_ > 0)

  // Step 1: y0 = L00^{-1} * x0
  Matrix x0(s0, nrhs);
  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < s0; ++i)
      x0(i, j) = x(i, j);
  p_child0_->applyInverseFactor(x0);

  if (isCholesky_ && rank_ > 0)
  {
    // Step 2: x1 = x1 - L10 * x0  where L10 = U_[1] * W_^T
    // temp = W_^T * x0  via dgemm  (rank_ x nrhs)
    int mT = static_cast<int>(rank_);
    int nT = static_cast<int>(nrhs);
    int kT = static_cast<int>(s0);
    double one = 1.0, zero = 0.0, neg = -1.0;
    int ldW = static_cast<int>(s0);
    int ldX0 = static_cast<int>(s0);
    int ldT = static_cast<int>(rank_);
    Matrix temp(mT, nT);
    HODLRDgemm("T", "N", &mT, &nT, &kT, &one, const_cast<double*>(&W_(0, 0)), &ldW, const_cast<double*>(&x0(0, 0)), &ldX0, &zero, &temp(0, 0), &ldT);

    // x1 -= U_[1] * temp  via dgemm
    int mU = static_cast<int>(s1);
    int nU = static_cast<int>(nrhs);
    int kU = static_cast<int>(rank_);
    int ldU1 = static_cast<int>(s1);
    int ldTemp = static_cast<int>(rank_);
    int ldX = static_cast<int>(size_);
    HODLRDgemm("N", "N", &mU, &nU, &kU, &neg, const_cast<double*>(&U_[1](0, 0)), &ldU1, &temp(0, 0), &ldTemp, &one, &x(s0, 0), &ldX);
  }

  // Step 3: x1 = L11^{-1} * x1  (recursive on child1)
  Matrix x1(s1, nrhs);
  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < s1; ++i)
      x1(i, j) = x(s0 + i, j);

  p_child1_->applyInverseFactor(x1);

  // Write back
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger i = 0; i < s0; ++i)
      x(i, j) = x0(i, j);
    for (UnsignedInteger i = 0; i < s1; ++i)
      x(s0 + i, j) = x1(i, j);
  }
}

void HODLRNode::applyInverseFactorTranspose(Matrix& x) const
{
  // Applies L^{-T} * x in-place where L is this node's HODLR Cholesky factor.
  // x must have size_ rows.
  const UnsignedInteger nrhs = x.getNbColumns();

  if (isLeaf_)
  {
    // Dense lower-triangular transpose solve via dtrsm: L^T * X = B
    int n = static_cast<int>(size_);
    int nrhs_int = static_cast<int>(nrhs);
    double one = 1.0;
    HODLRDtrsm("L", "L", "T", "N", &n, &nrhs_int, &one,
           const_cast<double*>(&Sfactor_(0, 0)), &n, &x(0, 0), &n);
    return;
  }

  const UnsignedInteger s0 = size_ / 2;
  const UnsignedInteger s1 = size_ - s0;

  // Backward substitution with the 2x2 block structure:
  // L^T = [L00^T, L10^T; 0, L11^T]
  // where L10 = U1 * W^T, so L10^T = W * U1^T
  //
  // 1. child1: L11^T * x1 = x1   (recursive)
  // 2. x0 = x0 - L10^T * x1 = x0 - W * U1^T * x1
  // 3. child0: L00^T * x0 = x0   (recursive)

  // Step 1: L11^T * x1 = x1  (recursive on child1)
  Matrix x1(s1, nrhs);
  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < s1; ++i)
      x1(i, j) = x(s0 + i, j);
  p_child1_->applyInverseFactorTranspose(x1);

  if (isCholesky_ && rank_ > 0)
  {
    // Step 2: x0 -= L10^T * x1  where L10^T = W * U1^T
    // temp = U1^T * x1  via dgemm  (rank_ x nrhs)
    int mU = static_cast<int>(rank_);
    int nU = static_cast<int>(nrhs);
    int kU = static_cast<int>(s1);
    double one = 1.0, zero = 0.0, neg = -1.0;
    int ldU1 = static_cast<int>(s1);
    int ldX1 = static_cast<int>(s1);
    int ldT = static_cast<int>(rank_);
    Matrix temp(mU, nU);
    HODLRDgemm("T", "N", &mU, &nU, &kU, &one, const_cast<double*>(&U_[1](0, 0)), &ldU1, const_cast<double*>(&x1(0, 0)), &ldX1, &zero, &temp(0, 0), &ldT);

    // x0 -= W * temp  via dgemm
    int mW = static_cast<int>(s0);
    int nW = static_cast<int>(nrhs);
    int kW = static_cast<int>(rank_);
    int ldW = static_cast<int>(s0);
    int ldTemp = static_cast<int>(rank_);
    int ldX0 = static_cast<int>(s0);
    HODLRDgemm("N", "N", &mW, &nW, &kW, &neg, const_cast<double*>(&W_(0, 0)), &ldW, &temp(0, 0), &ldTemp, &one, &x(0, 0), &ldX0);
  }

  // Step 3: L00^T * x0 = x0  (recursive on child0)
  Matrix x0(s0, nrhs);
  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < s0; ++i)
      x0(i, j) = x(i, j);
  p_child0_->applyInverseFactorTranspose(x0);

  // Write back
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger i = 0; i < s0; ++i)
      x(i, j) = x0(i, j);
    for (UnsignedInteger i = 0; i < s1; ++i)
      x(s0 + i, j) = x1(i, j);
  }
}

void HODLRNode::solve(Matrix& x) const
{
  HODLRBlasGuard blasGuard;
  if (isCholesky_)
  {
    // Cholesky: A = L * L^T, solve via forward + backward substitution
    applyInverseFactor(x);        // y = L^{-1} * x
    applyInverseFactorTranspose(x); // x = L^{-T} * y
    return;
  }

  // LU: recursive solve children then Woodbury correction
  if (!isLeaf_)
  {
    p_child0_->solve(x);
    p_child1_->solve(x);
  }
  applyInverse(x, 0);
}

void HODLRNode::solveLower(Matrix& x, Bool trans) const
{
  HODLRBlasGuard blasGuard;
  if (trans)
    applyInverseFactorTranspose(x);
  else
    applyInverseFactor(x);
}

Scalar HODLRNode::dotSolve(Matrix& x) const
{
  Matrix b(x);
  solve(b);
  Scalar result = 0.0;
  for (UnsignedInteger j = 0; j < x.getNbColumns(); ++j)
    for (UnsignedInteger i = 0; i < x.getNbRows(); ++i)
      result += x(i, j) * b(i, j);
  return result;
}

void HODLRNode::apply(Matrix& y, const Matrix& x) const
{
  HODLRBlasGuard blasGuard;
  const UnsignedInteger nrhs = x.getNbColumns();

  if (isLeaf_)
  {
    // y[start_:start_+size_] += leafMatrix_ * x[start_:start_+size_]
    const Scalar* leaf_data = &(*leafMatrix_.getImplementation())[0];
    for (UnsignedInteger j = 0; j < nrhs; ++j)
    {
      for (UnsignedInteger i = 0; i < size_; ++i)
      {
        Scalar sum = 0.0;
        for (UnsignedInteger k = 0; k < size_; ++k)
          sum += leaf_data[i + k * size_] * x(start_ + k, j);
        y(start_ + i, j) += sum;
      }
    }
    return;
  }

  const UnsignedInteger s0 = size_ / 2;
  const UnsignedInteger s1 = size_ - s0;

  // Recurse on children for diagonal blocks
  p_child0_->apply(y, x);
  p_child1_->apply(y, x);

  if (rank_ == 0)
    return;

  const Scalar* V0_data = &(*V_[0].getImplementation())[0];
  const Scalar* U1orig_data = &(*Uorig_[1].getImplementation())[0];

  // Off-diagonal: A[child1, child0] = Uorig_[1] * V_[0]^T
  // A[child0, child1] = V_[0] * Uorig_[1]^T  (symmetric)
  // temp0 = V_[0]^T * x[start_:start_+s0]     (rank_ x nrhs)
  // temp1 = Uorig_[1]^T * x[start_+s0:...]     (rank_ x nrhs)
  Matrix temp0(rank_, nrhs, 0.0);
  Matrix temp1(rank_, nrhs, 0.0);

  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger r = 0; r < rank_; ++r)
    {
      Scalar s0_val = 0.0;
      Scalar s1_val = 0.0;
      for (UnsignedInteger k = 0; k < s0; ++k)
        s0_val += V0_data[k + r * s0] * x(start_ + k, j);
      for (UnsignedInteger k = 0; k < s1; ++k)
        s1_val += U1orig_data[k + r * s1] * x(start_ + s0 + k, j);
      temp0(r, j) = s0_val;
      temp1(r, j) = s1_val;
    }
  }

  // y[start_+s0:...] += Uorig_[1] * temp0
  // y[start_:start_+s0] += V_[0] * temp1
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger r = 0; r < s1; ++r)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger k = 0; k < rank_; ++k)
        sum += U1orig_data[r + k * s1] * temp0(k, j);
      y(start_ + s0 + r, j) += sum;
    }
    for (UnsignedInteger r = 0; r < s0; ++r)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger k = 0; k < rank_; ++k)
        sum += V0_data[r + k * s0] * temp1(k, j);
      y(start_ + r, j) += sum;
    }
  }
}

void HODLRNode::applyFactor(Matrix& y, const Matrix& x) const
{
  HODLRBlasGuard blasGuard;
  const UnsignedInteger nrhs = x.getNbColumns();

  if (isLeaf_)
  {
    // y += L * x where L is the lower-triangular Cholesky factor in Sfactor_
    // For each RHS column: dtrmv does y += L * x
    for (UnsignedInteger j = 0; j < nrhs; ++j)
    {
      // Copy x vector for dtrmv (it writes result in-place)
      std::vector<double> tmp(size_);
      for (UnsignedInteger i = 0; i < size_; ++i)
        tmp[i] = x(start_ + i, j);
      int n = static_cast<int>(size_);
      int incx = 1;
      HODLRDtrmv("L", "N", "N", &n, const_cast<double*>(&Sfactor_(0, 0)), &n, tmp.data(), &incx);
      for (UnsignedInteger i = 0; i < size_; ++i)
        y(start_ + i, j) += tmp[i];
    }
    return;
  }

  const UnsignedInteger s0 = size_ / 2;
  const UnsignedInteger s1 = size_ - s0;

  // Recurse on children for diagonal blocks: L00, L11
  p_child0_->applyFactor(y, x);
  p_child1_->applyFactor(y, x);

  if (rank_ == 0 || !isCholesky_)
    return;

  const Scalar* W_data = &(*W_.getImplementation())[0];
  const Scalar* U1_data = &(*U_[1].getImplementation())[0];

  // Off-diagonal: L[child1, child0] = U_[1] * W_^T
  // y[start_+s0:...] += U_[1] * (W_^T * x[start_:start_+s0])
  // temp = W_^T * x0   (rank_ x nrhs)
  Matrix temp(rank_, nrhs, 0.0);
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger r = 0; r < rank_; ++r)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger k = 0; k < s0; ++k)
        sum += W_data[k + r * s0] * x(start_ + k, j);
      temp(r, j) = sum;
    }
  }

  // y[s0:s0+s1] += U_[1] * temp
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger r = 0; r < s1; ++r)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger k = 0; k < rank_; ++k)
        sum += U1_data[r + k * s1] * temp(k, j);
      y(start_ + s0 + r, j) += sum;
    }
  }
}

void HODLRNode::applyFactorTranspose(Matrix& y, const Matrix& x) const
{
  HODLRBlasGuard blasGuard;
  const UnsignedInteger nrhs = x.getNbColumns();

  if (isLeaf_)
  {
    // y += L^T * x where L is the lower-triangular Cholesky factor in Sfactor_
    for (UnsignedInteger j = 0; j < nrhs; ++j)
    {
      std::vector<double> tmp(size_);
      for (UnsignedInteger i = 0; i < size_; ++i)
        tmp[i] = x(start_ + i, j);
      int n = static_cast<int>(size_);
      int incx = 1;
      HODLRDtrmv("L", "T", "N", &n, const_cast<double*>(&Sfactor_(0, 0)), &n, tmp.data(), &incx);
      for (UnsignedInteger i = 0; i < size_; ++i)
        y(start_ + i, j) += tmp[i];
    }
    return;
  }

  const UnsignedInteger s0 = size_ / 2;
  const UnsignedInteger s1 = size_ - s0;

  // Recurse on children for diagonal blocks: L00^T, L11^T
  p_child0_->applyFactorTranspose(y, x);
  p_child1_->applyFactorTranspose(y, x);

  if (rank_ == 0 || !isCholesky_)
    return;

  const Scalar* W_data = &(*W_.getImplementation())[0];
  const Scalar* U1_data = &(*U_[1].getImplementation())[0];

  // Off-diagonal: L^T[child0, child1] = (U_[1] * W_^T)^T = W_ * U_[1]^T
  // y[start_:start_+s0] += W_ * (U_[1]^T * x[start_+s0:...])
  // temp = U_[1]^T * x1   (rank_ x nrhs)
  Matrix temp(rank_, nrhs, 0.0);
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger r = 0; r < rank_; ++r)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger k = 0; k < s1; ++k)
        sum += U1_data[k + r * s1] * x(start_ + s0 + k, j);
      temp(r, j) = sum;
    }
  }

  // y[0:s0] += W_ * temp
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger r = 0; r < s0; ++r)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger k = 0; k < rank_; ++k)
        sum += W_data[r + k * s0] * temp(k, j);
      y(start_ + r, j) += sum;
    }
  }
}

size_t HODLRNode::getNnz() const
{
  if (isLeaf_)
    return static_cast<size_t>(size_) * size_;
  else
  {
    size_t nnz = static_cast<size_t>(U_[0].getNbRows()) * U_[0].getNbColumns()
               + static_cast<size_t>(U_[1].getNbRows()) * U_[1].getNbColumns();
    if (p_child0_) nnz += p_child0_->getNnz();
    if (p_child1_) nnz += p_child1_->getNnz();
    return nnz;
  }
}

END_NAMESPACE_OPENTURNS
