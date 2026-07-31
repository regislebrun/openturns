%feature("docstring") OT::HODLRMatrix
"HODLR compressed matrix.

This class provides a compressed representation of dense matrices
using Hierarchical Off-Diagonal Low-Rank (HODLR) structure.
It is an alternative to :class:`~openturns.HMatrix` that does not require
external dependencies (hmat-oss or STRUMPACK).

Notes
-----
This class is experimental."

%feature("docstring") OT::HODLRMatrix::getNbRows
"Return the number of rows.

Returns
-------
n : int
    Number of rows."

%feature("docstring") OT::HODLRMatrix::getNbColumns
"Return the number of columns.

Returns
-------
n : int
    Number of columns."

%feature("docstring") OT::HODLRMatrix::factorize
"Factorize the matrix.

Parameters
----------
method : str
    Factorization method, one of 'LU', 'LLt' or 'LLT'."

%feature("docstring") OT::HODLRMatrix::setPermutation
"Set the spatial permutation.

The matrix is assembled in the permuted order, so this method must be
called before :meth:`assemble`. The permutation is a bijection of
``[0, ..., n-1]`` mapping the permuted index to the original one; an
empty permutation restores the original order.

Parameters
----------
permutation : sequence of int
    Permutation of ``[0, ..., n-1]``, or an empty sequence for the
    identity."

%feature("docstring") OT::HODLRMatrix::getPermutation
"Return the spatial permutation.

Returns
-------
permutation : :class:`~openturns.Indices`
    The permutation used at assembly time, empty for the identity."

%feature("docstring") OT::HODLRMatrix::solve
"Solve a linear system.

Parameters
----------
b : 1-d or 2-d sequence of float
    Right-hand side.
trans : bool, optional
    Whether to solve with the transpose. The transposed solve is not
    implemented: only False is accepted, and True raises an error.
    Default is False.

Returns
-------
x : 1-d or 2-d sequence of float
    Solution."

%feature("docstring") OT::HODLRMatrix::logDeterminant
"Return the log-determinant.

The matrix must be factorized first.

Returns
-------
logdet : float
    Logarithm of the absolute value of the determinant."

%feature("docstring") OT::HODLRMatrix::applyNugget
"Add a relative nugget to the diagonal.

Adds a small multiple of the mean diagonal to the matrix, following
the formulation C = sigma^2 I + K used in the HODLR Gaussian-process
literature. The relative factor is read from the resource key
'HODLRMatrix-Nugget' (default 1.0e-8); setting it to zero disables
the nugget. The nugget bounds the condition number of ill-conditioned
covariance matrices (e.g. long-correlation kernels) and is applied
through the same mechanism as :meth:`addIdentity`.

Notes
-----
This method must be called after assembly and before factorization."

%feature("docstring") OT::HODLRMatrix::copy
"Copy the matrix.

Returns
-------
m : :class:`~openturns.HODLRMatrix`
    A copy of the matrix with its own internal structure."

%feature("docstring") OT::HODLRMatrix::assemble
"Assemble the matrix from an evaluator.

Parameters
----------
f : callable
    Entry evaluator, such that f(i, j) returns the entry (i, j) of the matrix.
symmetry : str
    Symmetry flag, one of 'L' (lower part only, symmetric matrix), 'U'
    (upper part only, symmetric matrix), 'N' (no symmetry) or 'F'
    (full matrix)."

%feature("docstring") OT::HODLRMatrix::assembleReal
"Assemble the matrix from a Python callable.

Parameters
----------
callable : callable
    Python function that takes two integers i and j and returns the entry
    (i, j) of the matrix.
symmetry : str
    Symmetry flag, one of 'L' (lower part only, symmetric matrix), 'U'
    (upper part only, symmetric matrix), 'N' (no symmetry) or 'F'
    (full matrix)."

%feature("docstring") OT::HODLRMatrix::scale
"Scale the matrix by a scalar.

The matrix is multiplied by the scalar alpha, and the scaling is applied
to the internal representation without reassembling the matrix.

Parameters
----------
alpha : float
    Scaling factor."

%feature("docstring") OT::HODLRMatrix::gemv
"Compute a matrix-vector product.

Computes y = alpha * A * x + beta * y, where A is the matrix. Only the
non-transposed product is supported.

Parameters
----------
trans : str
    Transpose flag, only 'N' is supported.
alpha : float
    Scalar multiplier for A.
x : sequence of float
    Input vector.
beta : float
    Scalar multiplier for the output vector.
y : sequence of float
    Output vector, updated in place."

%feature("docstring") OT::HODLRMatrix::addIdentity
"Add a multiple of the identity to the matrix.

The matrix A is replaced by A + alpha * I. The matrix must be factorized
again after this call.

Parameters
----------
alpha : float
    Shift to add to the diagonal."

%feature("docstring") OT::HODLRMatrix::norm
"Return the matrix norm.

This method is not implemented yet."

%feature("docstring") OT::HODLRMatrix::getDiagonal
"Return the diagonal of the matrix.

Returns
-------
diagonal : :class:`~openturns.Point`
    Diagonal of the matrix, of dimension n."

%feature("docstring") OT::HODLRMatrix::solveLower
"Solve a linear system with the Cholesky factor.

Solves L * x = b (or L^T * x = b if trans is True), where L is the lower
Cholesky factor of the matrix. The matrix must have been factorized with
the 'LLt' or 'LLT' method.

Parameters
----------
b : 1-d or 2-d sequence of float
    Right-hand side.
trans : bool, optional
    Whether to solve with the transpose of the factor. Default is False.

Returns
-------
x : 1-d or 2-d sequence of float
    Solution."

%feature("docstring") OT::HODLRMatrix::compressionRatio
"Return the compression ratio of the matrix.

Returns
-------
ratio : sequence of int
    Pair (compressedSize, uncompressedSize), where compressedSize is the
    number of stored coefficients and uncompressedSize is the number of
    entries of the dense matrix."

%feature("docstring") OT::HODLRMatrix::fullrkRatio
"Return the full-rank ratio of the matrix.

Returns
-------
ratio : sequence of int
    Pair (fullRankSize, lowRankSize), where fullRankSize is the number of
    entries of the dense matrix and lowRankSize is the total storage of
    the low-rank blocks."

%feature("docstring") OT::HODLRMatrix::dump
"Write a debugging dump of the matrix.

Parameters
----------
name : str
    Name of the dump."
