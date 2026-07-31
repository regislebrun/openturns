//                                               -*- C++ -*-
/**
 *  @file  HODLRCore.hxx
 *  @brief Internal HODLR recursive tree using LAPACK
 *
 *  Based on the algorithm from george (Dan Foreman-Mackey, MIT license).
 *  Extended with OpenTURNS Matrix and LAPACK backend.
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
#ifndef OPENTURNS_HODLRCORE_HXX
#define OPENTURNS_HODLRCORE_HXX

#include "openturns/OTprivate.hxx"
#include "openturns/Matrix.hxx"
#include "openturns/Pointer.hxx"

#include <random>
#include <vector>
#include <functional>

BEGIN_NAMESPACE_OPENTURNS

/**
 * RAII guard that pins BLAS to single-thread mode.
 * HODLR operates on small matrices (rank=10, leaf=32-256) where
 * multi-threaded BLAS overhead dominates performance.
 * Restores the previous thread count on destruction.
 */
class OT_API HODLRBlasGuard
{
public:
  HODLRBlasGuard();
  ~HODLRBlasGuard();
private:
  int savedNumThreads_ = 1;
};

/**
 * Abstract interface for matrix entry evaluation.
 */
class OT_API HODLREntryEvaluator
{
public:
  virtual ~HODLREntryEvaluator() {}
  virtual HODLREntryEvaluator* clone() const
  {
    throw NotYetImplementedException(HERE) << "HODLREntryEvaluator::clone";
  }
  virtual Scalar operator()(UnsignedInteger, UnsignedInteger) const
  {
    throw NotYetImplementedException(HERE) << "HODLREntryEvaluator::operator()";
  }
  virtual UnsignedInteger getSize() const
  {
    throw NotYetImplementedException(HERE) << "HODLREntryEvaluator::getSize";
  }
};

/**
 * Evaluator that wraps an original evaluator and subtracts low-rank corrections.
 * Used for hierarchical Schur complement in HODLR Cholesky factorization.
 *
 * Stores the ROOT original evaluator and a flat list of all accumulated corrections
 * from every level of the recursion, avoiding nested evaluator chains that compound
 * numerical errors.
 *
 * corrected_eval(i, j) = root_original(i, j) - sum_k (U1_k * K_k * U1_k^T) + sum_k (lambda_k * delta(i, j))
 */
class OT_API HODLRCorrectedEvaluator : public HODLREntryEvaluator
{
public:
  /** Data for a single low-rank correction term */
  struct Correction {
    UnsignedInteger offset;
    UnsignedInteger size;
    UnsignedInteger rank;
    Matrix U1;   // n_rows x rank
    Matrix UK;   // U1 * K (precomputed, n_rows x rank)
    Scalar lambda;
  };

  /** Create a flattened evaluator with one or more corrections */
  HODLRCorrectedEvaluator(Pointer<const HODLREntryEvaluator> original,
                           UnsignedInteger offset,
                           UnsignedInteger size,
                           const std::vector<Correction>& corrections)
    : original_(original)
    , offset_(offset)
    , size_(size)
    , corrections_(corrections)
  {
  }

  /** Create a flattened evaluator from an existing one + a new correction */
  static Pointer<HODLRCorrectedEvaluator> flatten(
      const HODLRCorrectedEvaluator& existing,
      const Matrix& newU1,
      const Matrix& newK,
      Scalar newLambda,
      UnsignedInteger newOffset,
      UnsignedInteger newSize);

  /** Create a flattened evaluator with a single correction (first level) */
  static Pointer<HODLRCorrectedEvaluator> create(
      Pointer<const HODLREntryEvaluator> original,
      UnsignedInteger offset,
      UnsignedInteger size,
      const Matrix& U1,
      const Matrix& K,
      Scalar lambda);

  Scalar operator()(UnsignedInteger i, UnsignedInteger j) const override
  {
    Scalar val = (*original_)(i, j);
    for (const auto& corr : corrections_)
    {
      const SignedInteger local_i = static_cast<SignedInteger>(i) - static_cast<SignedInteger>(corr.offset);
      const SignedInteger local_j = static_cast<SignedInteger>(j) - static_cast<SignedInteger>(corr.offset);
      if (local_i >= 0 && local_j >= 0 &&
          static_cast<UnsignedInteger>(local_i) < corr.size &&
          static_cast<UnsignedInteger>(local_j) < corr.size)
      {
        Scalar term = 0.0;
        for (UnsignedInteger q = 0; q < corr.rank; ++q)
          term += corr.UK(local_i, q) * corr.U1(local_j, q);
        val -= term;
        if (i == j) val += corr.lambda;
      }
    }
    return val;
  }

  UnsignedInteger getSize() const override
  {
    return original_->getSize();
  }

  const Pointer<const HODLREntryEvaluator> getOriginal() const { return original_; }
  const std::vector<Correction>& getCorrections() const { return corrections_; }
  UnsignedInteger getOffset() const { return offset_; }
  UnsignedInteger getSizeLocal() const { return size_; }
  UnsignedInteger getNumCorrections() const { return corrections_.size(); }

private:
  Pointer<const HODLREntryEvaluator> original_;  // ALWAYS root original
  UnsignedInteger offset_;
  UnsignedInteger size_;
  std::vector<Correction> corrections_;

  // Friend declarations for internal access
  friend class HODLRNode;
};

/**
 * Internal HODLR tree node (symmetric matrices only).
 *
 * Implements the recursive HODLR structure from george (MIT license).
 * A01 = U * V^T, A10 = V * U^T.
 */
class OT_API HODLRNode
{
public:
  HODLRNode(Pointer<const HODLREntryEvaluator> eval,
            const Scalar* diag,
            UnsignedInteger start,
            UnsignedInteger size,
            UnsignedInteger minLeafSize,
            UnsignedInteger maxRank,
            Scalar tolerance,
            std::mt19937& rng,
            SignedInteger direction = 0,
            HODLRNode* parent = nullptr);

  ~HODLRNode();

  HODLRNode(const HODLRNode&) = delete;
  HODLRNode& operator=(const HODLRNode&) = delete;

  void compute();
  void computeCholesky();
  void solve(Matrix& x) const;
  Scalar dotSolve(Matrix& x) const;
  void solveLower(Matrix& x, Bool trans) const;  // L^{-1} or L^{-T} times x for the Cholesky factor L
  void apply(Matrix& y, const Matrix& x) const;
  void applyFactor(Matrix& y, const Matrix& x) const;

  Scalar getLogDeterminant() const { return logDet_; }
  UnsignedInteger getTotalRank() const { return totalRank_; }
  UnsignedInteger getNumLeaves() const { return numLeaves_; }
  size_t getNnz() const;
  void setShift(Scalar shift);

private:
  UnsignedInteger lowRankApprox(UnsignedInteger startRow, UnsignedInteger nRows,
                                UnsignedInteger startCol, UnsignedInteger nCols,
                                Scalar tol, std::mt19937& rng,
                                Matrix& Uout, Matrix& Vout);
  void factorize();
  void factorizeLeafCholesky();
  void factorizeLeafCholeskyCorrected(const Matrix& K, const Matrix& U1, Scalar lambda);  // Bypass evaluator, use dgemm correction
  void applyInverse(Matrix& x, UnsignedInteger start) const;
  void applyInverseLeaf(Matrix& x, UnsignedInteger start) const;
  void applyInverseCholeskyLeaf(Matrix& x, UnsignedInteger start) const;
  void applyInverseFactor(Matrix& x) const;  // L^{-1} * x for HODLR Cholesky factor L
  void applyInverseFactorTranspose(Matrix& x) const;  // L^{-T} * x for HODLR Cholesky factor L

  const Scalar* p_diag_;
  Pointer<const HODLREntryEvaluator> p_eval_;
  UnsignedInteger start_;
  UnsignedInteger size_;
  SignedInteger direction_;
  UnsignedInteger rank_;
  UnsignedInteger maxRank_;
  UnsignedInteger minLeafSize_;
  UnsignedInteger denseThreshold_;
  Scalar tolerance_;
  std::mt19937* p_rng_; // non-owning pointer to shared rng (evolves across the tree)
  bool isLeaf_;
  bool isCholesky_;
  Scalar logDet_;
  Scalar shift_;
  UnsignedInteger totalRank_;
  UnsignedInteger numLeaves_;

  HODLRNode* p_parent_;         // non-owning back-reference
  Pointer<HODLRNode> p_child0_;
  Pointer<HODLRNode> p_child1_;

  Collection<Matrix> U_;
  Collection<Matrix> V_;

  Matrix Sfactor_;
  Matrix leafMatrix_;
  Collection<Matrix> Uorig_;
  Matrix W_;   // L00^{-1} * V_[0] for internal nodes (Cholesky mode)
  // std::vector for LAPACK ABI compatibility (mutable int* required)
  std::vector<int> ipiv_;
};

END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_HODLRCORE_HXX */
