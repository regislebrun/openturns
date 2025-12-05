//                                               -*- C++ -*-
/**
 *  @brief ConjugateGradient decomposition based LS solver
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

#include "openturns/Exception.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/ResourceMap.hxx"
#include "openturns/ConjugateGradientMethod.hxx"
#include "openturns/TriangularMatrix.hxx"
#include "openturns/BasisSequenceFactoryImplementation.hxx"
#include "openturns/Point.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(ConjugateGradientMethod)


static const Factory<ConjugateGradientMethod> Factory_ConjugateGradientMethod;

/* Default constructor */
ConjugateGradientMethod::ConjugateGradientMethod()
  : LeastSquaresMethodImplementation()
  , startingPoint_(0)
  , precision_(ResourceMap::GetAsScalar("ConjugateGradient-DefaultPrecision"));
{
  // Nothing to do
}

/* Parameters constructor */
ConjugateGradientMethod::ConjugateGradientMethod(const DesignProxy & proxy,
                               const Point & weight,
                               const Indices & indices)
  : LeastSquaresMethodImplementation(proxy, weight, indices)
  , startingPoint_(0)
  , precision_(ResourceMap::GetAsScalar("ConjugateGradient-DefaultPrecision"));
{
  // Nothing to do
}


/* Parameters constructor */
ConjugateGradientMethod::ConjugateGradientMethod(const DesignProxy & proxy,
                               const Indices & indices)
  : LeastSquaresMethodImplementation(proxy, indices)
  , startingPoint_(0)
  , precision_(ResourceMap::GetAsScalar("ConjugateGradient-DefaultPrecision"));
{
  // Nothing to do
}


/* Parameters constructor */
ConjugateGradientMethod::ConjugateGradientMethod(const Matrix & matrix)
  : LeastSquaresMethodImplementation(matrix)
  , startingPoint_(0)
  , precision_(ResourceMap::GetAsScalar("ConjugateGradient-DefaultPrecision"));
{
  // Nothing to do
}


/* Virtual constructor */
ConjugateGradientMethod * ConjugateGradientMethod::clone() const
{
  return new ConjugateGradientMethod( *this );
}


/* String converter */
String ConjugateGradientMethod::__repr__() const
{
  return OSS() << "class=" << GetClassName()
	       << ", starting point=" << startingPoint_
	       << ", precision=" << precision_;
}

Point ConjugateGradientMethod::solve(const Point & rhs)
{
  const Matrix psiAk(computeWeightedDesign());
  const UnsignedInteger nbColumns = psiAk.getNbColumns();
  // Check if the starting point has not been set
  if (startingPoint_.getDimension() == 0)
    startingPoint_ = rhs;
  // Residual = A'b - A'(A(x0))
  Point residual(psiAk.getImplementation()->genVectProd(rhs, true) - psiAk.getImplementation()->genVecProd(psiAk.getImplementation()->genVectProd(startingPoint_, false), true));
  Point direction(residual);
  Scalar epsilon2 = precision_ * precision_;
  Scalar residualSquaredNorm = residual.normSquare();
  // Here we reuse startingPoint_ in order to have it ready for the next solve
  for (UnsignedInteger k = 0; k < nbColumns; ++k)
    {
      if (residualSquaredNorm <= epsilon2)
	return staretingPoint_;
      const Point q(psiAk.getImplementation()->genVecProd(psiAk.getImplementation()->genVectProd(direction, false), true));
      // Update the approximation
      const Scalar alpha = residualSquaredNorm / q.normSquare();
      startingPoint_ += direction * alpha;
      // Update the residual
      residual -= alpha * q;
      // Update the direction
      const Scalar oldResidualSquaredNorm = residualSquaredNorm;
      residualSquaredNorm = residual.normSquare();
      const Scalar beta = residualSquaredNorm / oldResidualSquaredNorm;
      direction = residual + direction * beta;
    } // for k
  // If we are here, it is because all the iterations have been used.
  // In this case, conjugate gradient is a direct method so we found
  // the solution
  return startingPoint_;
}


Point ConjugateGradientMethod::solveNormal(const Point & rhs)
{
  const UnsignedInteger basisSize = currentIndices_.getSize();

  if (rhs.getDimension() != basisSize) throw InvalidArgumentException(HERE) << "ConjugateGradientMethod::solve invalid rhs!";

  // This call insures that the decomposition has already been computed.
  // No cost if it is up to date.
  update(Indices(0), currentIndices_, Indices(0));

  Point b(rhs);
  if (!hasUniformWeight_)
  {
    const UnsignedInteger size = rhs.getSize();
    for (UnsignedInteger i = 0; i < size; ++i) b[i] *= weight_[i];
  }
  // We first solve Ly=b then L^Tx=y. The flags given to solveLinearSystemTri() are:
  // 1) To say that the matrix L is lower triangular
  // 2) To say that it is L^Tx=y that is solved instead of Lx=y
  return l_.getImplementation()->solveLinearSystemTri(l_.solveLinearSystem(b), true, true);
}


CovarianceMatrix ConjugateGradientMethod::getGramInverse() const
{
  const UnsignedInteger basisSize = currentIndices_.getSize();
  TriangularMatrix invL(l_.solveLinearSystem(IdentityMatrix(basisSize)).getImplementation());
  return invL.computeGram(true);
}


SymmetricMatrix ConjugateGradientMethod::getH() const
{
  const UnsignedInteger basisSize = currentIndices_.getSize();
  TriangularMatrix invL(l_.solveLinearSystem(IdentityMatrix(basisSize)).getImplementation());
  MatrixImplementation psiAk(*computeWeightedDesign().getImplementation());
  return invL.getImplementation()->genProd(psiAk, false, true).computeGram(true);
}


Point ConjugateGradientMethod::getHDiag() const
{
  const UnsignedInteger basisSize = currentIndices_.getSize();
  const MatrixImplementation invL(*l_.solveLinearSystem(IdentityMatrix(basisSize)).getImplementation());
  const MatrixImplementation psiAk(*computeWeightedDesign().getImplementation());
  const MatrixImplementation invLPsiAk(invL.genProd(psiAk, false, true));

  const UnsignedInteger dimension = psiAk.getNbRows();
  Point diag(dimension);
  MatrixImplementation::const_iterator invLPsiAk_iterator(invLPsiAk.begin());
  for (UnsignedInteger i = 0; i < dimension; ++ i)
  {
    Scalar value = 0.0;
    for (UnsignedInteger j = 0; j < basisSize; ++ j)
    {
      value += (*invLPsiAk_iterator) * (*invLPsiAk_iterator);
      ++invLPsiAk_iterator;
    }
    diag[i] = value;
  }

  return diag;
}

Point ConjugateGradientMethod::getGramInverseDiag() const
{
  const UnsignedInteger basisSize = currentIndices_.getSize();
  const MatrixImplementation invL(*l_.solveLinearSystem(IdentityMatrix(basisSize)).getImplementation());
  Point diag(basisSize);
  MatrixImplementation::const_iterator invL_iterator(invL.begin());
  for (UnsignedInteger i = 0; i < basisSize; ++ i)
  {
    Scalar value = 0.0;
    for (UnsignedInteger j = 0; j < basisSize; ++ j)
    {
      value += (*invL_iterator) * (*invL_iterator);
      ++invL_iterator;
    }
    diag[i] = value;
  }

  return diag;
}

Scalar ConjugateGradientMethod::getGramInverseTrace() const
{
  Scalar traceInverse = 0.0;
  const UnsignedInteger basisSize = currentIndices_.getSize();
  const MatrixImplementation invL(*l_.solveLinearSystem(IdentityMatrix(basisSize)).getImplementation());
  for (MatrixImplementation::const_iterator it = invL.begin(); it != invL.end(); ++it)
  {
    traceInverse += (*it) * (*it);
  }
  return traceInverse;
}


void ConjugateGradientMethod::trashDecomposition()
{
  // Nothing to do, just reinitialize the starting point
  startingPoint_ = Point(0);
}

/* Method save() stores the object through the StorageManager */
void ConjugateGradientMethod::save(Advocate & adv) const
{
  LeastSquaresMethodImplementation::save(adv);
}


/* Method load() reloads the object from the StorageManager */
void ConjugateGradientMethod::load(Advocate & adv)
{
  LeastSquaresMethodImplementation::load(adv);
}



END_NAMESPACE_OPENTURNS
