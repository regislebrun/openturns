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
#include <cmath>
#include <chrono>
#include <cstdio>

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

// Temporary factorization timing instrumentation (removed after tuning).
struct FactorTiming
{
  double recomp = 0.0;
  double recompCopy = 0.0;
  double recompQR = 0.0;
  double recompMcat = 0.0;
  double recompGram = 0.0;
  double recompSyevd = 0.0;
  double recompFinal = 0.0;
  double leaf = 0.0;
  double leafCorr = 0.0;
  double leafCorrected = 0.0;
  double W = 0.0;
  double Kgram = 0.0;
  double UK = 0.0;
  double dpotrf = 0.0;
  int nRecomp = 0;
  int nLeaf = 0;
  int nLeafCorr = 0;
  int nLeafCorrected = 0;
  double k1sum = 0.0;
  double nRowsSum = 0.0;
  double totalRankSum = 0.0;
  ~FactorTiming()
  {
    if (printProfile)
      std::fprintf(stderr, "[FT] recomp: %.2f ms (%d) [copy %.2f QR %.2f Mcat %.2f Gram %.2f syevd %.2f final %.2f; k1 %.0f nRows %.0f totRank %.0f], leaf: %.2f ms (%d) of which corrections %.2f ms, leafCorrected: %.2f ms (%d), W: %.2f, Kgram: %.2f, UK: %.2f, dpotrf: %.2f\n",
                   (recompCopy + recompQR + recompMcat + recompGram + recompSyevd + recompFinal) * 1e3, nRecomp,
                   recompCopy * 1e3, recompQR * 1e3, recompMcat * 1e3,
                   recompGram * 1e3, recompSyevd * 1e3, recompFinal * 1e3,
                   k1sum, nRowsSum, totalRankSum,
                   leaf * 1e3, nLeaf, leafCorr * 1e3, leafCorrected * 1e3, nLeafCorrected,
                   W * 1e3, Kgram * 1e3, UK * 1e3, dpotrf * 1e3);
  }
  // Set once during factorization: only print the profile when explicitly
  // requested, to avoid polluting stderr in normal use.
  bool printProfile = false;
} g_factorTiming;

}  // namespace

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

// Fill a symmetric dense leaf block using only the lower-triangular kernel
// evaluations and mirror them into the upper triangle. The HODLR tree is
// symmetric by construction, and the entry evaluators used (kernel,
// permuted, corrected) all satisfy f(i,j) = f(j,i), so the mirrored values
// are exact.
void HODLRSymmetricLeafFill(MatrixImplementation& Sfact,
                            UnsignedInteger size,
                            UnsignedInteger start,
                            const HODLREntryEvaluator& eval)
{
  for (UnsignedInteger j = 0; j < size; ++j)
  {
    for (UnsignedInteger i = j; i < size; ++i)
    {
      const Scalar v = eval(start + i, start + j);
      Sfact[i + j * size] = v;
      if (i != j)
        Sfact[j + i * size] = v;
    }
  }
}

}  // namespace

HODLRNode::HODLRNode(Pointer<const HODLREntryEvaluator> eval,
                     const Scalar* diag,
                     UnsignedInteger start,
                     UnsignedInteger size,
                     UnsignedInteger minLeafSize,
                     UnsignedInteger maxRank,
                     Scalar tolerance,
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
    // Max-element (partial-pivoting) ACA.
    rank_ = lowRankApproxPartialPivot(
        start_ + half, size_ - half,
        start_, half,
        tolerance, U_[1], V_[0]);

    // Fall back to dense when the block is not low-rank: if the ACA reached
    // the full rank min(rows, cols), store the block exactly as U_[1] = A01
    // and V_[0] = I instead of keeping a truncated (and inexact) product.
    const UnsignedInteger s1 = size_ - half;
    if (rank_ == std::min(s1, half))
    {
      U_[1] = Matrix(s1, half);
      V_[0] = Matrix(half, half);
      MatrixImplementation& U1Impl = *U_[1].getImplementation();
      MatrixImplementation& V0Impl = *V_[0].getImplementation();
      for (UnsignedInteger j = 0; j < half; ++j)
        V0Impl[j + j * half] = 1.0;
      for (UnsignedInteger j = 0; j < half; ++j)
        for (UnsignedInteger i = 0; i < s1; ++i)
          U1Impl[i + j * s1] = (*p_eval_)(start_ + half + i, start_ + j);
    }

    // Symmetric: A10 = A01^T = V_[0] * U_[1]^T, so:
    U_[0] = V_[0];
    V_[1] = U_[1];

    totalRank_ = rank_;

    p_child0_ = new HODLRNode(p_eval_, p_diag_, start_, half, minLeafSize, maxRank_, tolerance_, 0, this);
    p_child1_ = new HODLRNode(p_eval_, p_diag_, start_ + half, size_ - half, minLeafSize, maxRank_, tolerance_, 1, this);
  }
  else
  {
    isLeaf_ = true;
    // Cache the raw kernel leaf block at construction time so that the
    // factorization (LU or Cholesky) never re-evaluates the kernel entries.
    leafKernel_ = Matrix(size_, size_);
    MatrixImplementation& leafData = *leafKernel_.getImplementation();
    HODLRSymmetricLeafFill(leafData, size_, start_, *p_eval_);
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

UnsignedInteger HODLRNode::lowRankApproxPartialPivot(UnsignedInteger startRow, UnsignedInteger nRows,
    UnsignedInteger startCol, UnsignedInteger nCols,
    Scalar tol, Matrix& Uout, Matrix& Vout)
{
  const UnsignedInteger maxRankBound = std::min(nRows, nCols);
  const UnsignedInteger maxRank = (maxRank_ > 0) ? std::min(maxRankBound, maxRank_) : maxRankBound;

  // The factor buffers are fully written below before any read (each new U/V
  // column fills every row), so the zero-initialization would only add a
  // useless memset of the whole block. They are kept in transposed layout
  // (factor row l contiguous in the block row/column index) so that the
  // residual downdates are contiguous rank-1 loops that the compiler can
  // vectorize, instead of dot products with a strided access.
  //
  // The buffers grow geometrically with the rank revealed by the ACA. The
  // worst-case allocation min(nRows, nCols) * (nRows + nCols) would otherwise
  // be paid on every block even when the rank stays small (e.g. ~100 MB for a
  // 2500 x 2500 block of an exactly low-rank 1D kernel), and the value
  // initialization of std::vector would zero it all.
  UnsignedInteger bufferCols = 0;
  std::vector<Scalar> Udata;
  std::vector<Scalar> Vdata;
  const auto ensureBuffer = [&](const UnsignedInteger needed)
  {
    if (needed <= bufferCols)
      return;
    const UnsignedInteger grown =
        std::max(needed, std::max(UnsignedInteger(1), 2 * bufferCols));
    Udata.resize(grown * nRows);
    Vdata.resize(grown * nCols);
    bufferCols = grown;
  };

  // Partial-pivoting ACA, ported from the "ACA partial" scheme of hmat-oss
  // (src/compression.cpp, doCompressionAcaPartial). Unlike the full-pivoting
  // variant, the residual block is never materialized: only the residual row
  // at the pivot row and the residual column at the pivot column are assembled
  // from the evaluator on demand, so the cost is O(rank * (nRows + nCols))
  // instead of O(rank * nRows * nCols). Two ingredients make the pivot choice
  // as reliable as full pivoting:
  //  - rows and columns already used as pivots are excluded (rowFree/colFree),
  //    which prevents the stagnation of naive partial pivoting, and
  //  - the pivot column is scanned over the residual of the pivot row only,
  //    while the next pivot row is the row with the largest residual in the
  //    pivot column.
  // The stopping criterion accumulates the same Frobenius-norm estimate as
  // the full-pivoting variant (without the cross terms between the factors).
  UnsignedInteger rank = 0;
  Scalar norm = 0.0;
  const Scalar tol2 = tol * tol;

  std::vector<unsigned char> rowFree(nRows, 1);
  std::vector<unsigned char> colFree(nCols, 1);
  std::vector<Scalar> aCol(nRows, 0.0);
  std::vector<Scalar> bCol(nCols, 0.0);
  UnsignedInteger pivotRow = 0;
  UnsignedInteger pivotCol = 0;

  while (rank < maxRank)
  {
    // Residual row at the pivot row: A(pivotRow, :) - sum_l U[pivotRow, l] * V[:, l]
    for (UnsignedInteger j = 0; j < nCols; ++j)
      bCol[j] = (*p_eval_)(startRow + pivotRow, startCol + j);
    for (UnsignedInteger l = 0; l < rank; ++l)
    {
      const Scalar uval = Udata[l * nRows + pivotRow];
      Scalar* const b = bCol.data();
      const Scalar* const vptr = Vdata.data() + l * nCols;
#pragma omp simd
      for (UnsignedInteger j = 0; j < nCols; ++j)
        b[j] -= uval * vptr[j];
    }
    Scalar maxVal = 0.0;
    for (UnsignedInteger j = 0; j < nCols; ++j)
    {
      const Scalar absVal = std::abs(bCol[j]);
      if (colFree[j] && (absVal > maxVal))
      {
        maxVal = absVal;
        pivotCol = j;
      }
    }
    rowFree[pivotRow] = 0;

    if (maxVal < 1e-14)
    {
      // The residual row is (numerically) zero on all free columns. Because
      // pivotRow is the free row with the largest residual in the previous
      // pivot column (see below), a zero residual there means the block is
      // exhausted: for a low-rank block the remaining rows are near-proportional
      // to it and would also be below the threshold. Terminate immediately
      // instead of re-evaluating every remaining row (which is O(nRows*nCols)
      // kernel evaluations on an already-captured block, e.g. an exactly
      // low-rank 1D block such as a degree-2 exponential-polynomial kernel).
      break;
    }

    // Residual column at the pivot column: A(:, pivotCol) - sum_l U[:, l] * V[pivotCol, l]
    for (UnsignedInteger i = 0; i < nRows; ++i)
      aCol[i] = (*p_eval_)(startRow + i, startCol + pivotCol);
    for (UnsignedInteger l = 0; l < rank; ++l)
    {
      const Scalar vval = Vdata[l * nCols + pivotCol];
      Scalar* const a = aCol.data();
      const Scalar* const uptr = Udata.data() + l * nRows;
#pragma omp simd
      for (UnsignedInteger i = 0; i < nRows; ++i)
        a[i] -= vval * uptr[i];
    }
    colFree[pivotCol] = 0;

    // Same scaling convention as the full-pivoting variant: the pivot column
    // is stored unscaled in U, the pivot row divided by the pivot in V.
    const Scalar pivot = bCol[pivotCol];
    ensureBuffer(rank + 1);
    Scalar uNorm2 = 0.0;
    {
      Scalar* const uptr = Udata.data() + rank * nRows;
      for (UnsignedInteger i = 0; i < nRows; ++i)
      {
        const Scalar v = aCol[i];
        uptr[i] = v;
        uNorm2 += v * v;
      }
    }
    Scalar vNorm2 = 0.0;
    {
      Scalar* const vptr = Vdata.data() + rank * nCols;
      const Scalar invPivot = 1.0 / pivot;
      for (UnsignedInteger j = 0; j < nCols; ++j)
      {
        const Scalar v = bCol[j] * invPivot;
        vptr[j] = v;
        vNorm2 += v * v;
      }
    }
    const Scalar rowcolNorm = uNorm2 * vNorm2;

    norm += rowcolNorm;
    ++rank;

    if (rowcolNorm < tol2 * norm) break;

    if (rank == maxRank) break;

    // Next pivot row: the row with the largest residual in the pivot column,
    // among the rows not yet used as pivots.
    maxVal = 0.0;
    pivotRow = 0;
    for (UnsignedInteger i = 0; i < nRows; ++i)
    {
      const Scalar absVal = std::abs(aCol[i]);
      if (rowFree[i] && (absVal > maxVal))
      {
        maxVal = absVal;
        pivotRow = i;
      }
    }
    if (maxVal < 1e-14)
      break;
  }

  if (rank == 0)
  {
    // The block is numerically zero: no pivot survives the hard-coded
    // threshold. Keep the deterministic rank-maxRank representation used
    // by lowRankApprox (U = first columns of A, V = I) so downstream code
    // never has to deal with an empty U/V pair.
    if ((maxRank_ > 0) && (maxRank_ < maxRankBound))
      LOGWARN(OSS() << "HODLRMatrix: rank-starved block (" << nRows << "x" << nCols
             << ") reached the max rank " << maxRank_ << " before the assembly tolerance; "
             << "accuracy may be degraded. Set HODLRMatrix-MaxRank to 0 for adaptive (tolerance-driven) rank");
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

  Uout = Matrix(nRows, rank);
  Vout = Matrix(nCols, rank);
  {
    MatrixImplementation& UoutImpl = *Uout.getImplementation();
    MatrixImplementation& VoutImpl = *Vout.getImplementation();
    Scalar* const Uout_data = &UoutImpl[0];
    Scalar* const Vout_data = &VoutImpl[0];
    for (UnsignedInteger l = 0; l < rank; ++l)
    {
      const Scalar* const up = Udata.data() + l * nRows;
      const Scalar* const vp = Vdata.data() + l * nCols;
      for (UnsignedInteger i = 0; i < nRows; ++i)
        Uout_data[i + l * nRows] = up[i];
      for (UnsignedInteger j = 0; j < nCols; ++j)
        Vout_data[j + l * nCols] = vp[j];
    }
  }

  if ((maxRank_ > 0) && (rank >= maxRank_) && (rank < maxRankBound))
    LOGWARN(OSS() << "HODLRMatrix: rank-starved block (" << nRows << "x" << nCols
           << ") reached the max rank " << maxRank_ << " before the assembly tolerance; "
           << "accuracy may be degraded. Set HODLRMatrix-MaxRank to 0 for adaptive (tolerance-driven) rank");

  return rank;
}

void HODLRNode::recompressLowRank(Matrix& Uout, Matrix& Vout,
    UnsignedInteger startRow, UnsignedInteger nRows,
    UnsignedInteger startCol, UnsignedInteger nCols,
    const std::vector<HODLRCorrectedEvaluator::Correction>& corrections,
    Scalar tolerance)
{
  // Recompress the low-rank block A01' = Uout * Vout^T - sum_k (UK_k * U1_k^T)
  // where Uout/Vout approximate the pure kernel block (from the initial ACA)
  // and each correction term is an exact low-rank matrix of the hierarchical
  // Schur complement. Instead of re-running ACA with a corrected evaluator,
  // the correction is subtracted algebraically and the sum of low-rank terms
  // is truncated back to the assembly tolerance by SVD.
  //
  // The block row range [startRow, startRow+nRows) and column range
  // [startCol, startCol+nCols) are assumed fully covered by every correction,
  // which is guaranteed by the Schur complement recursion of computeCholesky.
  const UnsignedInteger rankIn = Uout.getNbColumns();
  if (rankIn == 0) return;

  const auto ft0 = std::chrono::steady_clock::now();
  ++g_factorTiming.nRecomp;
  g_factorTiming.nRowsSum += nRows;

  UnsignedInteger totalRank = rankIn;
  for (const auto& corr : corrections)
    totalRank += corr.rank;
  g_factorTiming.totalRankSum += totalRank;

  // Ucat: nRows x totalRank, Vcat: nCols x totalRank
  MatrixImplementation Ucat(nRows, totalRank);
  MatrixImplementation Vcat(nCols, totalRank);
  {
    const Scalar* Uout_data = &(*Uout.getImplementation())[0];
    const Scalar* Vout_data = &(*Vout.getImplementation())[0];
    const UnsignedInteger ldUout = Uout.getNbRows();
    const UnsignedInteger ldVout = Vout.getNbRows();
    Scalar* Ucat_data = &Ucat[0];
    Scalar* Vcat_data = &Vcat[0];
    for (UnsignedInteger c = 0; c < rankIn; ++c)
    {
      for (UnsignedInteger r = 0; r < nRows; ++r)
        Ucat_data[r + c * nRows] = Uout_data[r + c * ldUout];
      for (UnsignedInteger r = 0; r < nCols; ++r)
        Vcat_data[r + c * nCols] = Vout_data[r + c * ldVout];
    }
    UnsignedInteger col = rankIn;
    for (const auto& corr : corrections)
    {
      if (corr.rank == 0) continue;
      const UnsignedInteger rowOffset = startRow - corr.offset;
      const UnsignedInteger colOffset = startCol - corr.offset;
      const Scalar* uk_data = &(*corr.UK.getImplementation())[0];
      const Scalar* u1_data = &(*corr.U1.getImplementation())[0];
      const UnsignedInteger ldCorr = corr.UK.getNbRows();
      for (UnsignedInteger c = 0; c < corr.rank; ++c)
      {
        for (UnsignedInteger r = 0; r < nRows; ++r)
          Ucat_data[r + col * nRows] = -uk_data[(rowOffset + r) + c * ldCorr];
        for (UnsignedInteger r = 0; r < nCols; ++r)
          Vcat_data[r + col * nCols] = u1_data[(colOffset + r) + c * ldCorr];
        ++col;
      }
    }
    totalRank = col;
  }
  g_factorTiming.recompCopy += std::chrono::duration<double>(std::chrono::steady_clock::now() - ft0).count();
  g_factorTiming.k1sum += std::min(nRows, totalRank);

  // QR of Ucat in place: Ucat = Q * R (thin QR, k1 = min(nRows, totalRank)).
  // Qcat overwrites Ucat; R is extracted as the leading k1 x totalRank upper triangle.
  const UnsignedInteger k1 = std::min(nRows, totalRank);
  MatrixImplementation Rcat(k1, totalRank);
  {
    const auto t0 = std::chrono::steady_clock::now();
    int m = static_cast<int>(nRows);
    int n = static_cast<int>(totalRank);
    int lda = static_cast<int>(nRows);
    int kq = static_cast<int>(k1);
    Point tau(k1);
    double lworkQuery = 0.0;
    int lwork = -1;
    int info = 0;
    dgeqrf_(&m, &n, &Ucat[0], &lda, &tau[0], &lworkQuery, &lwork, &info);
    lwork = static_cast<int>(lworkQuery);
    Point work(lwork);
    dgeqrf_(&m, &n, &Ucat[0], &lda, &tau[0], &work[0], &lwork, &info);
    if (info != 0)
      throw InternalException(HERE) << "dgeqrf failed (info=" << info << ")";
    // Extract R (upper trapezoid) BEFORE dorgqr overwrites Ucat with Q.
    {
      const Scalar* Ud = &Ucat[0];
      Scalar* Rd = &Rcat[0];
      for (UnsignedInteger c = 0; c < totalRank; ++c)
      {
        const UnsignedInteger lim = std::min(c + 1, k1);
        for (UnsignedInteger r = 0; r < lim; ++r)
          Rd[r + c * k1] = Ud[r + c * nRows];
      }
    }
    lwork = -1;
    dorgqr_(&m, &kq, &kq, &Ucat[0], &lda, &tau[0], &lworkQuery, &lwork, &info);
    lwork = static_cast<int>(lworkQuery);
    work = Point(lwork);
    dorgqr_(&m, &kq, &kq, &Ucat[0], &lda, &tau[0], &work[0], &lwork, &info);
    if (info != 0)
      throw InternalException(HERE) << "dorgqr failed (info=" << info << ")";
    g_factorTiming.recompQR += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }
  // SVD of M = Vcat * R^T (nCols x k1) through its Gram matrix G = M^T * M.
  // Two algebraically equivalent routes:
  //  - direct:  form M, then G = M^T M              (cost O(nCols*k1*totalRank))
  //  - merged:  H = Vcat^T Vcat, G = R * H * R^T    (cost O(nCols*totalRank^2))
  // The merged route avoids the large nCols x k1 intermediate and is cheaper
  // when the block is much wider than the rank (nCols > 2*k1).
  const bool useMergedGram = (nCols > 2 * k1);
  MatrixImplementation G(k1, k1);
  MatrixImplementation Mcat;  // only allocated in the direct route
  {
    const auto t0 = std::chrono::steady_clock::now();
    if (useMergedGram)
    {
      MatrixImplementation H(totalRank, totalRank);
      {
        int mm = static_cast<int>(totalRank);
        int nn = static_cast<int>(totalRank);
        int kk = static_cast<int>(nCols);
        double one = 1.0;
        double zero = 0.0;
        int ldV = static_cast<int>(nCols);
        int ldH = static_cast<int>(totalRank);
        HODLRDgemm("T", "N", &mm, &nn, &kk, &one, &Vcat[0], &ldV, &Vcat[0], &ldV, &zero, &H[0], &ldH);
      }
      MatrixImplementation Tmp(totalRank, k1);
      {
        int mm = static_cast<int>(totalRank);
        int nn = static_cast<int>(k1);
        int kk = static_cast<int>(totalRank);
        double one = 1.0;
        double zero = 0.0;
        int ldH = static_cast<int>(totalRank);
        int ldR = static_cast<int>(k1);
        int ldT = static_cast<int>(totalRank);
        HODLRDgemm("N", "T", &mm, &nn, &kk, &one, &H[0], &ldH, &Rcat[0], &ldR, &zero, &Tmp[0], &ldT);
      }
      {
        int mm = static_cast<int>(k1);
        int nn = static_cast<int>(k1);
        int kk = static_cast<int>(totalRank);
        double one = 1.0;
        double zero = 0.0;
        int ldR = static_cast<int>(k1);
        int ldT = static_cast<int>(totalRank);
        int ldG = static_cast<int>(k1);
        HODLRDgemm("N", "N", &mm, &nn, &kk, &one, &Rcat[0], &ldR, &Tmp[0], &ldT, &zero, &G[0], &ldG);
      }
    }
    else
    {
      Mcat = MatrixImplementation(nCols, k1);
      {
        int mm = static_cast<int>(nCols);
        int nn = static_cast<int>(k1);
        int kk = static_cast<int>(totalRank);
        double one = 1.0;
        double zero = 0.0;
        int ldV = static_cast<int>(nCols);
        int ldR = static_cast<int>(k1);
        int ldM = static_cast<int>(nCols);
        HODLRDgemm("N", "T", &mm, &nn, &kk, &one, &Vcat[0], &ldV, &Rcat[0], &ldR, &zero, &Mcat[0], &ldM);
      }
      {
        int mm = static_cast<int>(k1);
        int nn = static_cast<int>(k1);
        int kk = static_cast<int>(nCols);
        double one = 1.0;
        double zero = 0.0;
        int ldM = static_cast<int>(nCols);
        int ldG = static_cast<int>(k1);
        HODLRDgemm("T", "N", &mm, &nn, &kk, &one, &Mcat[0], &ldM, &Mcat[0], &ldM, &zero, &G[0], &ldG);
      }
    }
    g_factorTiming.recompGram += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }

  Point eigenvalues(k1);
  {
    const auto t0 = std::chrono::steady_clock::now();
    int n = static_cast<int>(k1);
    char jobz = 'V';
    char uplo = 'L';
    int ljobz = 1;
    int luplo = 1;
    int lwork = -1;
    int liwork = -1;
    double lworkQuery = 0.0;
    int liworkQuery = 0;
    int info = 0;
    dsyevd_(&jobz, &uplo, &n, &G[0], &n, &eigenvalues[0], &lworkQuery, &lwork, &liworkQuery, &liwork, &info, &ljobz, &luplo);
    lwork = static_cast<int>(lworkQuery);
    liwork = static_cast<int>(liworkQuery);
    Point work(lwork);
    std::vector<int> iwork(static_cast<std::size_t>(liwork));
    dsyevd_(&jobz, &uplo, &n, &G[0], &n, &eigenvalues[0], &work[0], &lwork, &iwork[0], &liwork, &info, &ljobz, &luplo);
    if (info != 0)
      throw InternalException(HERE) << "dsyevd failed to converge (info=" << info << ")";
    g_factorTiming.recompSyevd += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }
  const UnsignedInteger nSing = std::min(k1, nCols);
  const UnsignedInteger maxRank = std::min(nRows, nCols);
  UnsignedInteger rankOut = std::min(nSing, maxRank);
  Scalar totalNorm2 = 0.0;
  for (UnsignedInteger i = 0; i < k1; ++i)
    if (eigenvalues[i] > 0.0)
      totalNorm2 += eigenvalues[i];
  if ((rankOut > 0) && (tolerance > 0.0) && (totalNorm2 > 0.0))
  {
    // Keep the largest r singular values such that the Frobenius norm of the
    // truncated part is below tolerance relative to the block norm (same
    // criterion as the ACA stopping rule).
    const Scalar tol2 = tolerance * tolerance;
    Scalar droppedNorm2 = 0.0;
    rankOut = nSing;
    for (UnsignedInteger i = k1 - nSing; i < k1; ++i)
    {
      const Scalar lambda = (eigenvalues[i] > 0.0) ? eigenvalues[i] : 0.0;
      droppedNorm2 += lambda;
      if (droppedNorm2 <= tol2 * totalNorm2)
        rankOut = k1 - 1 - i;
      else
        break;
    }
    if (rankOut == 0) rankOut = 1;
    if (rankOut > maxRank) rankOut = maxRank;
  }

  // Eigenvectors live in the columns of G (ascending eigenvalues); take the
  // rankOut largest, in descending order (columns k1-1 down to k1-rankOut),
  // scaled by s and 1/s so that U_new = Qcat * Ps, V_new = Mcat * Pi (direct
  // route) or V_new = Vcat * (Rcat^T * Pi) (merged route).
  MatrixImplementation Unew(nRows, rankOut);
  MatrixImplementation Vnew(nCols, rankOut);
  if (rankOut > 0)
  {
    const auto t0 = std::chrono::steady_clock::now();
    const Scalar* G_data = &G[0];
    MatrixImplementation Ps(k1, rankOut);
    MatrixImplementation Pi(k1, rankOut);
    Scalar* Ps_data = &Ps[0];
    Scalar* Pi_data = &Pi[0];
    for (UnsignedInteger c = 0; c < rankOut; ++c)
    {
      const UnsignedInteger col = k1 - 1 - c;
      const Scalar s = std::sqrt(eigenvalues[col] > 0.0 ? eigenvalues[col] : 0.0);
      const Scalar invS = (s > 0.0) ? (1.0 / s) : 0.0;
      for (UnsignedInteger r = 0; r < k1; ++r)
      {
        const Scalar p = G_data[r + col * k1];
        Ps_data[r + c * k1] = p * s;
        Pi_data[r + c * k1] = p * invS;
      }
    }
    {
      int mm = static_cast<int>(nRows);
      int nn = static_cast<int>(rankOut);
      int kk = static_cast<int>(k1);
      double one = 1.0;
      double zero = 0.0;
      int ldQ = static_cast<int>(nRows);
      int ldP = static_cast<int>(k1);
      int ldU = static_cast<int>(nRows);
      HODLRDgemm("N", "N", &mm, &nn, &kk, &one, &Ucat[0], &ldQ, &Ps[0], &ldP, &zero, &Unew[0], &ldU);
    }
    if (useMergedGram)
    {
      MatrixImplementation RPi(totalRank, rankOut);
      {
        int mm = static_cast<int>(totalRank);
        int nn = static_cast<int>(rankOut);
        int kk = static_cast<int>(k1);
        double one = 1.0;
        double zero = 0.0;
        int ldR = static_cast<int>(k1);
        int ldP = static_cast<int>(k1);
        int ldT = static_cast<int>(totalRank);
        HODLRDgemm("T", "N", &mm, &nn, &kk, &one, &Rcat[0], &ldR, &Pi[0], &ldP, &zero, &RPi[0], &ldT);
      }
      {
        int mm = static_cast<int>(nCols);
        int nn = static_cast<int>(rankOut);
        int kk = static_cast<int>(totalRank);
        double one = 1.0;
        double zero = 0.0;
        int ldV = static_cast<int>(nCols);
        int ldT = static_cast<int>(totalRank);
        int ldVn = static_cast<int>(nCols);
        HODLRDgemm("N", "N", &mm, &nn, &kk, &one, &Vcat[0], &ldV, &RPi[0], &ldT, &zero, &Vnew[0], &ldVn);
      }
    }
    else
    {
      int mm = static_cast<int>(nCols);
      int nn = static_cast<int>(rankOut);
      int kk = static_cast<int>(k1);
      double one = 1.0;
      double zero = 0.0;
      int ldM = static_cast<int>(nCols);
      int ldP = static_cast<int>(k1);
      int ldV = static_cast<int>(nCols);
      HODLRDgemm("N", "N", &mm, &nn, &kk, &one, &Mcat[0], &ldM, &Pi[0], &ldP, &zero, &Vnew[0], &ldV);
    }
    g_factorTiming.recompFinal += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }

  Uout = Matrix(nRows, rankOut);
  Vout = Matrix(nCols, rankOut);
  {
    MatrixImplementation& UoutImpl = *Uout.getImplementation();
    MatrixImplementation& VoutImpl = *Vout.getImplementation();
    std::copy(&Unew[0], &Unew[0] + nRows * rankOut, &UoutImpl[0]);
    std::copy(&Vnew[0], &Vnew[0] + nCols * rankOut, &VoutImpl[0]);
  }
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
  g_factorTiming.printProfile = ResourceMap::GetAsBool("HODLRMatrix-ProfileFactorization");
  computeCholesky(std::vector<HODLRCorrectedEvaluator::Correction>());
}

void HODLRNode::computeCholesky(const std::vector<HODLRCorrectedEvaluator::Correction>& corrections)
{
  HODLRBlasGuard blasGuard;
  isCholesky_ = true;

  if (isLeaf_)
  {
    const auto t0 = std::chrono::steady_clock::now();
    factorizeLeafCholesky(corrections);
    g_factorTiming.leaf += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    ++g_factorTiming.nLeaf;
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

  const UnsignedInteger s0 = size_ / 2;
  const UnsignedInteger s1 = size_ - s0;

  // Correct this node's off-diagonal block A01' = A01 - sum corrections.
  // The kernel approximation factors computed at assembly time are reused and
  // recompressed by SVD, avoiding a full re-assembly with a corrected evaluator.
  if (!corrections.empty() && ResourceMap::GetAsBool("HODLRMatrix-RecompressCorrections"))
  {
    const auto t0 = std::chrono::steady_clock::now();
    recompressLowRank(U_[1], V_[0], start_ + s0, s1, start_, s0, corrections, tolerance_);
    g_factorTiming.recomp += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    rank_ = U_[1].getNbColumns();
    U_[0] = V_[0];
    V_[1] = U_[1];
  }

  // 1. Recursively compute L00 = chol(A00')
  p_child0_->computeCholesky(corrections);

  totalRank_ = rank_ + p_child0_->getTotalRank();
  numLeaves_ = p_child0_->getNumLeaves();

  if (rank_ > 0)
  {
    // 2. Compute W = L00^{-1} * V_[0]
    W_ = V_[0];
    {
      const auto t0 = std::chrono::steady_clock::now();
      p_child0_->applyInverseFactor(W_);
      g_factorTiming.W += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }

    // 3. Build K = W^T * W  (rank_ x rank_) via dgemm
    Matrix K(rank_, rank_);
    {
      int mK = static_cast<int>(rank_);
      int nK = static_cast<int>(rank_);
      int kK = static_cast<int>(s0);
      double one = 1.0, zero = 0.0;
      int ldW = static_cast<int>(s0);
      int ldK = static_cast<int>(rank_);
      const auto t0 = std::chrono::steady_clock::now();
      HODLRDgemm("T", "N", &mK, &nK, &kK, &one, &W_(0, 0), &ldW, &W_(0, 0), &ldW, &zero, &K(0, 0), &ldK);
      g_factorTiming.Kgram += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
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
                                  s1 + 1, maxRank_, tolerance_, 1, this);
        p_child1_->setShift(shift_);
        try
        {
          p_child1_->factorizeLeafCholeskyCorrected(K, U_[1], lambda, corrections);
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
        // Hierarchical Schur complement: propagate the new correction to child1's
        // subtree in place. The factors assembled with the plain kernel evaluator
        // are recompressed by SVD at each level (see recompressLowRank), so no
        // tree rebuild with a corrected evaluator is required.
        std::vector<HODLRCorrectedEvaluator::Correction> childCorrections = corrections;
        HODLRCorrectedEvaluator::Correction newCorr;
        newCorr.offset = start_ + s0;
        newCorr.size = s1;
        newCorr.rank = rank_;
        newCorr.U1 = U_[1];
        {
          Matrix UK(s1, rank_);
          int m = static_cast<int>(s1);
          int k = static_cast<int>(rank_);
          int l = static_cast<int>(rank_);
          double one = 1.0, zero = 0.0;
          const auto t0 = std::chrono::steady_clock::now();
          HODLRDgemm("N", "N", &m, &k, &l, &one,
                 const_cast<double*>(&U_[1](0, 0)), &m,
                 const_cast<double*>(&K(0, 0)), &l,
                 &zero, &UK(0, 0), &m);
          g_factorTiming.UK += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
          newCorr.UK = UK;
        }
        newCorr.lambda = lambda;
        childCorrections.push_back(newCorr);
        try
        {
          p_child1_->computeCholesky(childCorrections);
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
    p_child1_->computeCholesky(corrections);
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
  g_factorTiming.printProfile = ResourceMap::GetAsBool("HODLRMatrix-ProfileFactorization");
  if (isLeaf_)
  {
    // Use the kernel block cached at construction; only the shift changes.
    // Deep-copy: Matrix assignment would alias leafKernel_ through the
    // shared implementation, and the in-place LU below must not destroy it.
    Sfactor_ = Matrix(*leafKernel_.getImplementation());
    {
      MatrixImplementation& Sfact = *Sfactor_.getImplementation();
      if (shift_ != 0.0)
        for (UnsignedInteger i = 0; i < size_; ++i)
          Sfact[i + i * size_] += shift_;
    }
    leafMatrix_ = Matrix(*Sfactor_.getImplementation());
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

void HODLRNode::factorizeLeafCholesky(const std::vector<HODLRCorrectedEvaluator::Correction>& corrections)
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

  // Leaf: assemble the corrected dense matrix from the kernel block cached
  // at construction (deep-copied so the in-place dgemm corrections and dpotrf
  // below do not destroy leafKernel_), then apply every accumulated Schur
  // complement correction via batch dgemm.
  Sfactor_ = Matrix(*leafKernel_.getImplementation());
  {
    MatrixImplementation& Sfact = *Sfactor_.getImplementation();

    // Apply global shift to all diagonals (from addIdentity)
    if (shift_ != 0.0)
      for (UnsignedInteger i = 0; i < size_; ++i)
        Sfact[i + i * size_] += shift_;

    // Apply each correction via batch dgemm
    const UnsignedInteger n = size_;
    const UnsignedInteger leafStart = start_;
    const auto t0 = std::chrono::steady_clock::now();
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
    g_factorTiming.leafCorr += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    ++g_factorTiming.nLeafCorr;
  }

  int info = 0;
  int n = static_cast<int>(size_);
  {
    MatrixImplementation& Sfact2 = *Sfactor_.getImplementation();
    const auto t0 = std::chrono::steady_clock::now();
    HODLRDpotrf("L", &n, &Sfact2[0], &n, &info);
    g_factorTiming.dpotrf += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }
  if (info != 0)
    throw InternalException(HERE) << "Cholesky factorization failed, info=" << info;
}

void HODLRNode::factorizeLeafCholeskyCorrected(const Matrix& K, const Matrix& U1, Scalar lambda,
                                               const std::vector<HODLRCorrectedEvaluator::Correction>& corrections)
{
  // Directly assemble A11' = A11 - sum corrections - U1 * K * U1^T + lambda * I
  // Uses the ORIGINAL evaluator (p_eval_) for A11 entries, computes
  // the low-rank corrections via dgemm for numerical stability.
  const UnsignedInteger n = size_;
  const UnsignedInteger rank = K.getNbRows();
  const auto ft0 = std::chrono::steady_clock::now();
  ++g_factorTiming.nLeafCorrected;

  Sfactor_ = Matrix(*leafKernel_.getImplementation());
  MatrixImplementation& Sfact = *Sfactor_.getImplementation();

  // 1b. Apply the accumulated Schur complement corrections
  const UnsignedInteger leafStart = start_;
  for (const auto& corr : corrections)
  {
    const SignedInteger localOffset =
        static_cast<SignedInteger>(corr.offset) - static_cast<SignedInteger>(leafStart);
    const SignedInteger localEnd =
        localOffset + static_cast<SignedInteger>(corr.size);
    if (localEnd <= 0) continue;
    if (localOffset >= static_cast<SignedInteger>(n)) break;

    const UnsignedInteger i0 = static_cast<UnsignedInteger>(std::max(SignedInteger(0), localOffset));
    const UnsignedInteger i1 = static_cast<UnsignedInteger>(std::min(static_cast<SignedInteger>(n), localEnd));
    const UnsignedInteger localN = i1 - i0;

    if (localN == 0 || corr.rank == 0) continue;

    const UnsignedInteger ukRowOffset = i0 - static_cast<UnsignedInteger>(localOffset);
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

    if (corr.lambda != 0.0)
      for (UnsignedInteger ii = i0; ii < i1; ++ii)
        Sfact[ii + ii * n] += corr.lambda;
  }

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
  g_factorTiming.leafCorrected += std::chrono::duration<double>(std::chrono::steady_clock::now() - ft0).count();
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
  const UnsignedInteger nxRows = x.getNbRows();
  x.copyOnWrite();
  Scalar* x_data = &(*x.getImplementation())[0];

  // temp[0:rank_, :] = V1^T * x[s0:s0+s1, :]
  // temp[rank_:2*rank_, :] = V0^T * x[0:s0, :]
  const UnsignedInteger ntRows = 2 * rank_;
  Matrix temp(ntRows, nrhs, 0.0);
  MatrixImplementation& tempImpl = *temp.getImplementation();
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger i = 0; i < rank_; ++i)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger r = 0; r < s1; ++r)
        sum += V1_data[r + i * s1] * x_data[(offset + s0 + r) + j * nxRows];
      tempImpl[i + j * ntRows] = sum;
    }
    for (UnsignedInteger i = 0; i < rank_; ++i)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger r = 0; r < s0; ++r)
        sum += V0_data[r + i * s0] * x_data[(offset + r) + j * nxRows];
      tempImpl[(rank_ + i) + j * ntRows] = sum;
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
    HODLRDgetrs(&trans, &n, &nr, &(*Scopy.getImplementation())[0], &n, ipivCopy.data(), &tempImpl[0], &ldb, &info);
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
        sum += U0_data[r + i * s0] * tempImpl[i + j * ntRows];
      x_data[(offset + r) + j * nxRows] -= sum;
    }
    for (UnsignedInteger r = 0; r < s1; ++r)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger i = 0; i < rank_; ++i)
        sum += U1_data[r + i * s1] * tempImpl[(rank_ + i) + j * ntRows];
      x_data[(offset + s0 + r) + j * nxRows] -= sum;
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

  const UnsignedInteger nxRows = x.getNbRows();
  x.copyOnWrite();
  Scalar* x_data = &(*x.getImplementation())[0];

  Matrix block(size_, nrhs);
  MatrixImplementation& blockImpl = *block.getImplementation();
  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < size_; ++i)
      blockImpl[i + j * size_] = x_data[(offset + i) + j * nxRows];

  Matrix Dcopy(Sfactor_);
  std::vector<int> ipivCopy(ipiv_.begin(), ipiv_.end());
  HODLRDgetrs(&trans, &n, &nr, &(*Dcopy.getImplementation())[0], &n, ipivCopy.data(), &blockImpl[0], &ldb, &info);
  if (info != 0)
    throw InternalException(HERE) << "Solve failed, info=" << info;

  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < size_; ++i)
      x_data[(offset + i) + j * nxRows] = blockImpl[i + j * size_];
}

void HODLRNode::applyInverseCholeskyLeaf(Matrix& x, UnsignedInteger start) const
{
  const UnsignedInteger nrhs = x.getNbColumns();
  const SignedInteger offset = static_cast<SignedInteger>(start_) - static_cast<SignedInteger>(start);
  int n = static_cast<int>(size_);
  int nr = static_cast<int>(nrhs);
  int ldb = n;
  int info = 0;

  const UnsignedInteger nxRows = x.getNbRows();
  x.copyOnWrite();
  Scalar* x_data = &(*x.getImplementation())[0];

  Matrix block(size_, nrhs);
  MatrixImplementation& blockImpl = *block.getImplementation();
  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < size_; ++i)
      blockImpl[i + j * size_] = x_data[(offset + i) + j * nxRows];

  // Solve L * L^T * block = block  =>  block = C^{-1} * block
  Matrix Lcopy(Sfactor_);
  HODLRDpotrs("L", &n, &nr, &(*Lcopy.getImplementation())[0], &n, &blockImpl[0], &ldb, &info);
  if (info != 0)
    throw InternalException(HERE) << "Solve failed, info=" << info;

  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < size_; ++i)
      x_data[(offset + i) + j * nxRows] = blockImpl[i + j * size_];
}

void HODLRNode::applyInverseFactor(Matrix& x) const
{
  // Applies L^{-1} * x in-place where L is this node's HODLR Cholesky factor.
  // x must have size_ rows.
  const UnsignedInteger nrhs = x.getNbColumns();
  const UnsignedInteger nxRows = x.getNbRows();
  x.copyOnWrite();
  Scalar* x_data = &(*x.getImplementation())[0];

  if (isLeaf_)
  {
    // Dense lower-triangular solve via dtrsm
    int n = static_cast<int>(size_);
    int nrhs_int = static_cast<int>(nrhs);
    double one = 1.0;
    HODLRDtrsm("L", "L", "N", "N", &n, &nrhs_int, &one,
           const_cast<double*>(&(*Sfactor_.getImplementation())[0]), &n, x_data, &n);
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
  MatrixImplementation& x0Impl = *x0.getImplementation();
  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < s0; ++i)
      x0Impl[i + j * s0] = x_data[i + j * nxRows];
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
    MatrixImplementation& tempImpl = *temp.getImplementation();
    HODLRDgemm("T", "N", &mT, &nT, &kT, &one, const_cast<double*>(&(*W_.getImplementation())[0]), &ldW, &x0Impl[0], &ldX0, &zero, &tempImpl[0], &ldT);

    // x1 -= U_[1] * temp  via dgemm
    int mU = static_cast<int>(s1);
    int nU = static_cast<int>(nrhs);
    int kU = static_cast<int>(rank_);
    int ldU1 = static_cast<int>(s1);
    int ldTemp = static_cast<int>(rank_);
    int ldX = static_cast<int>(size_);
    HODLRDgemm("N", "N", &mU, &nU, &kU, &neg, const_cast<double*>(&(*U_[1].getImplementation())[0]), &ldU1, &tempImpl[0], &ldTemp, &one, &x_data[s0], &ldX);
  }

  // Step 3: x1 = L11^{-1} * x1  (recursive on child1)
  Matrix x1(s1, nrhs);
  MatrixImplementation& x1Impl = *x1.getImplementation();
  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < s1; ++i)
      x1Impl[i + j * s1] = x_data[(s0 + i) + j * nxRows];

  p_child1_->applyInverseFactor(x1);

  // Write back
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger i = 0; i < s0; ++i)
      x_data[i + j * nxRows] = x0Impl[i + j * s0];
    for (UnsignedInteger i = 0; i < s1; ++i)
      x_data[(s0 + i) + j * nxRows] = x1Impl[i + j * s1];
  }
}

void HODLRNode::applyInverseFactorTranspose(Matrix& x) const
{
  // Applies L^{-T} * x in-place where L is this node's HODLR Cholesky factor.
  // x must have size_ rows.
  const UnsignedInteger nrhs = x.getNbColumns();
  const UnsignedInteger nxRows = x.getNbRows();
  x.copyOnWrite();
  Scalar* x_data = &(*x.getImplementation())[0];

  if (isLeaf_)
  {
    // Dense lower-triangular transpose solve via dtrsm: L^T * X = B
    int n = static_cast<int>(size_);
    int nrhs_int = static_cast<int>(nrhs);
    double one = 1.0;
    HODLRDtrsm("L", "L", "T", "N", &n, &nrhs_int, &one,
           const_cast<double*>(&(*Sfactor_.getImplementation())[0]), &n, x_data, &n);
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
  MatrixImplementation& x1Impl = *x1.getImplementation();
  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < s1; ++i)
      x1Impl[i + j * s1] = x_data[(s0 + i) + j * nxRows];
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
    MatrixImplementation& tempImpl = *temp.getImplementation();
    HODLRDgemm("T", "N", &mU, &nU, &kU, &one, const_cast<double*>(&(*U_[1].getImplementation())[0]), &ldU1, &x1Impl[0], &ldX1, &zero, &tempImpl[0], &ldT);

    // x0 -= W * temp  via dgemm
    int mW = static_cast<int>(s0);
    int nW = static_cast<int>(nrhs);
    int kW = static_cast<int>(rank_);
    int ldW = static_cast<int>(s0);
    int ldTemp = static_cast<int>(rank_);
    int ldX0 = static_cast<int>(s0);
    HODLRDgemm("N", "N", &mW, &nW, &kW, &neg, const_cast<double*>(&(*W_.getImplementation())[0]), &ldW, &tempImpl[0], &ldTemp, &one, x_data, &ldX0);
  }

  // Step 3: L00^T * x0 = x0  (recursive on child0)
  Matrix x0(s0, nrhs);
  MatrixImplementation& x0Impl = *x0.getImplementation();
  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < s0; ++i)
      x0Impl[i + j * s0] = x_data[i + j * nxRows];
  p_child0_->applyInverseFactorTranspose(x0);

  // Write back
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger i = 0; i < s0; ++i)
      x_data[i + j * nxRows] = x0Impl[i + j * s0];
    for (UnsignedInteger i = 0; i < s1; ++i)
      x_data[(s0 + i) + j * nxRows] = x1Impl[i + j * s1];
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
  const UnsignedInteger nrhs = x.getNbColumns();
  const UnsignedInteger nxRows = x.getNbRows();
  const Scalar* x_data = &(*x.getImplementation())[0];
  const Scalar* b_data = &(*b.getImplementation())[0];
  for (UnsignedInteger j = 0; j < nrhs; ++j)
    for (UnsignedInteger i = 0; i < nxRows; ++i)
      result += x_data[i + j * nxRows] * b_data[i + j * nxRows];
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
    const UnsignedInteger nxRows = x.getNbRows();
    const Scalar* x_data = &(*x.getImplementation())[0];
    y.copyOnWrite();
    const UnsignedInteger nyRows = y.getNbRows();
    Scalar* y_data = &(*y.getImplementation())[0];
    for (UnsignedInteger j = 0; j < nrhs; ++j)
    {
      for (UnsignedInteger i = 0; i < size_; ++i)
      {
        Scalar sum = 0.0;
        for (UnsignedInteger k = 0; k < size_; ++k)
          sum += leaf_data[i + k * size_] * x_data[(start_ + k) + j * nxRows];
        y_data[(start_ + i) + j * nyRows] += sum;
      }
    }
    return;
  }

  const UnsignedInteger s0 = size_ / 2;
  const UnsignedInteger s1 = size_ - s0;
  const Scalar* V0_data = &(*V_[0].getImplementation())[0];
  const Scalar* U1orig_data = &(*Uorig_[1].getImplementation())[0];
  const UnsignedInteger nxRows = x.getNbRows();
  const Scalar* x_data = &(*x.getImplementation())[0];
  y.copyOnWrite();
  const UnsignedInteger nyRows = y.getNbRows();
  Scalar* y_data = &(*y.getImplementation())[0];

  // Recurse on children for diagonal blocks
  p_child0_->apply(y, x);
  p_child1_->apply(y, x);

  if (rank_ == 0)
    return;

  // Off-diagonal: A[child1, child0] = Uorig_[1] * V_[0]^T
  // A[child0, child1] = V_[0] * Uorig_[1]^T  (symmetric)
  // temp0 = V_[0]^T * x[start_:start_+s0]     (rank_ x nrhs)
  // temp1 = Uorig_[1]^T * x[start_+s0:...]     (rank_ x nrhs)
  Matrix temp0(rank_, nrhs, 0.0);
  Matrix temp1(rank_, nrhs, 0.0);
  MatrixImplementation& temp0Impl = *temp0.getImplementation();
  MatrixImplementation& temp1Impl = *temp1.getImplementation();

  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger r = 0; r < rank_; ++r)
    {
      Scalar s0_val = 0.0;
      Scalar s1_val = 0.0;
      for (UnsignedInteger k = 0; k < s0; ++k)
        s0_val += V0_data[k + r * s0] * x_data[(start_ + k) + j * nxRows];
      for (UnsignedInteger k = 0; k < s1; ++k)
        s1_val += U1orig_data[k + r * s1] * x_data[(start_ + s0 + k) + j * nxRows];
      temp0Impl[r + j * rank_] = s0_val;
      temp1Impl[r + j * rank_] = s1_val;
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
        sum += U1orig_data[r + k * s1] * temp0Impl[k + j * rank_];
      y_data[(start_ + s0 + r) + j * nyRows] += sum;
    }
    for (UnsignedInteger r = 0; r < s0; ++r)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger k = 0; k < rank_; ++k)
        sum += V0_data[r + k * s0] * temp1Impl[k + j * rank_];
      y_data[(start_ + r) + j * nyRows] += sum;
    }
  }
}

void HODLRNode::applyFactor(Matrix& y, const Matrix& x) const
{
  HODLRBlasGuard blasGuard;
  const UnsignedInteger nrhs = x.getNbColumns();
  const UnsignedInteger nxRows = x.getNbRows();
  const Scalar* x_data = &(*x.getImplementation())[0];
  y.copyOnWrite();
  const UnsignedInteger nyRows = y.getNbRows();
  Scalar* y_data = &(*y.getImplementation())[0];

  if (isLeaf_)
  {
    // y += L * x where L is the lower-triangular Cholesky factor in Sfactor_
    // For each RHS column: dtrmv does y += L * x
    for (UnsignedInteger j = 0; j < nrhs; ++j)
    {
      // Copy x vector for dtrmv (it writes result in-place)
      std::vector<double> tmp(size_);
      for (UnsignedInteger i = 0; i < size_; ++i)
        tmp[i] = x_data[(start_ + i) + j * nxRows];
      int n = static_cast<int>(size_);
      int incx = 1;
      HODLRDtrmv("L", "N", "N", &n, const_cast<double*>(&(*Sfactor_.getImplementation())[0]), &n, tmp.data(), &incx);
      for (UnsignedInteger i = 0; i < size_; ++i)
        y_data[(start_ + i) + j * nyRows] += tmp[i];
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
  MatrixImplementation& tempImpl = *temp.getImplementation();
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger r = 0; r < rank_; ++r)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger k = 0; k < s0; ++k)
        sum += W_data[k + r * s0] * x_data[(start_ + k) + j * nxRows];
      tempImpl[r + j * rank_] = sum;
    }
  }

  // y[s0:s0+s1] += U_[1] * temp
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger r = 0; r < s1; ++r)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger k = 0; k < rank_; ++k)
        sum += U1_data[r + k * s1] * tempImpl[k + j * rank_];
      y_data[(start_ + s0 + r) + j * nyRows] += sum;
    }
  }
}

void HODLRNode::applyFactorTranspose(Matrix& y, const Matrix& x) const
{
  HODLRBlasGuard blasGuard;
  const UnsignedInteger nrhs = x.getNbColumns();
  const UnsignedInteger nxRows = x.getNbRows();
  const Scalar* x_data = &(*x.getImplementation())[0];
  y.copyOnWrite();
  const UnsignedInteger nyRows = y.getNbRows();
  Scalar* y_data = &(*y.getImplementation())[0];

  if (isLeaf_)
  {
    // y += L^T * x where L is the lower-triangular Cholesky factor in Sfactor_
    for (UnsignedInteger j = 0; j < nrhs; ++j)
    {
      std::vector<double> tmp(size_);
      for (UnsignedInteger i = 0; i < size_; ++i)
        tmp[i] = x_data[(start_ + i) + j * nxRows];
      int n = static_cast<int>(size_);
      int incx = 1;
      HODLRDtrmv("L", "T", "N", &n, const_cast<double*>(&(*Sfactor_.getImplementation())[0]), &n, tmp.data(), &incx);
      for (UnsignedInteger i = 0; i < size_; ++i)
        y_data[(start_ + i) + j * nyRows] += tmp[i];
    }
    return;
  }

  const UnsignedInteger s0 = size_ / 2;
  const UnsignedInteger s1 = size_ - s0;

  // Recurse on children for the diagonal blocks: L00^T, L11^T
  p_child0_->applyFactorTranspose(y, x);
  p_child1_->applyFactorTranspose(y, x);

  if (rank_ == 0 || !isCholesky_)
    return;

  const Scalar* W_data = &(*W_.getImplementation())[0];
  const Scalar* U1_data = &(*U_[1].getImplementation())[0];

  // Off-diagonal: L[child1, child0] = U_[1] * W_^T, so the transposed block
  // is L[child0, child1]^T = W_ * U_[1]^T:
  // y[start_:start_+s0] += W_ * (U_[1]^T * x[start_+s0:...])
  // temp = U_[1]^T * x1   (rank_ x nrhs)
  Matrix temp(rank_, nrhs, 0.0);
  MatrixImplementation& tempImpl = *temp.getImplementation();
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger r = 0; r < rank_; ++r)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger k = 0; k < s1; ++k)
        sum += U1_data[k + r * s1] * x_data[(start_ + s0 + k) + j * nxRows];
      tempImpl[r + j * rank_] = sum;
    }
  }

  // y[s0:s0+s1]... no: y[0:s0] += W_ * temp
  for (UnsignedInteger j = 0; j < nrhs; ++j)
  {
    for (UnsignedInteger r = 0; r < s0; ++r)
    {
      Scalar sum = 0.0;
      for (UnsignedInteger k = 0; k < rank_; ++k)
        sum += W_data[r + k * s0] * tempImpl[k + j * rank_];
      y_data[(start_ + r) + j * nyRows] += sum;
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
