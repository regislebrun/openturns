%define OT_HODLRMatrixImplementation_doc
"HODLR compressed matrix implementation.

This class implements the compressed HODLR (Hierarchical Off-Diagonal Low-Rank)
matrix representation based on the algorithm from george (MIT license).
It uses randomized ACA for low-rank approximation and recursive LDLT/LU factorization.

Notes
-----
This class is experimental."
%enddef

%feature("docstring") OT::HODLRMatrixImplementation
OT_HODLRMatrixImplementation_doc

%feature("docstring") OT::HODLRMatrixImplementation::getNbRows
"Return the number of rows.

Returns
-------
n : int
    Number of rows."

%feature("docstring") OT::HODLRMatrixImplementation::getNbColumns
"Return the number of columns.

Returns
-------
n : int
    Number of columns."

%feature("docstring") OT::HODLRMatrixImplementation::getParameters
"Return the HODLR matrix parameters.

Returns
-------
parameters : :class:`~openturns.HODLRMatrixParameters`
    Parameters of the HODLR matrix."

%feature("docstring") OT::HODLRMatrixImplementation::assemble
"Assemble the matrix from an evaluator.

Parameters
----------
f : callable
    Entry evaluator, such that f(i, j) returns the entry (i, j) of the matrix.
symmetry : str
    Symmetry flag, one of 'L' (lower part only, symmetric matrix), 'U'
    (upper part only, symmetric matrix), 'N' (no symmetry) or 'F'
    (full matrix)."

%feature("docstring") OT::HODLRMatrixImplementation::factorize
"Factorize the matrix.

Parameters
----------
method : str
    Factorization method, one of 'LU', 'LLt' or 'LLT'."

%feature("docstring") OT::HODLRMatrixImplementation::scale
"Scale the matrix by a scalar.

The matrix is multiplied by the scalar alpha, and the scaling is applied
to the internal representation without reassembling the matrix.

Parameters
----------
alpha : float
    Scaling factor."

%feature("docstring") OT::HODLRMatrixImplementation::gemv
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

%feature("docstring") OT::HODLRMatrixImplementation::addIdentity
"Add a multiple of the identity to the matrix.

The matrix A is replaced by A + alpha * I. The matrix must be factorized
again after this call.

Parameters
----------
alpha : float
    Shift to add to the diagonal."

%feature("docstring") OT::HODLRMatrixImplementation::applyNugget
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

%feature("docstring") OT::HODLRMatrixImplementation::norm
"Return the matrix norm.

This method is not implemented yet."

%feature("docstring") OT::HODLRMatrixImplementation::getDiagonal
"Return the diagonal of the matrix.

Returns
-------
diagonal : :class:`~openturns.Point`
    Diagonal of the matrix, of dimension n."

%feature("docstring") OT::HODLRMatrixImplementation::solve
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

%feature("docstring") OT::HODLRMatrixImplementation::solveLower
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

%feature("docstring") OT::HODLRMatrixImplementation::logDeterminant
"Return the log-determinant.

The matrix must be factorized first.

Returns
-------
logdet : float
    Logarithm of the absolute value of the determinant."

%feature("docstring") OT::HODLRMatrixImplementation::compressionRatio
"Return the compression ratio of the matrix.

Returns
-------
ratio : sequence of int
    Pair (compressedSize, uncompressedSize), where compressedSize is the
    number of stored coefficients and uncompressedSize is the number of
    entries of the dense matrix."

%feature("docstring") OT::HODLRMatrixImplementation::fullrkRatio
"Return the full-rank ratio of the matrix.

Returns
-------
ratio : sequence of int
    Pair (fullRankSize, lowRankSize), where fullRankSize is the number of
    entries of the dense matrix and lowRankSize is the total storage of
    the low-rank blocks."

%feature("docstring") OT::HODLRMatrixImplementation::dump
"Write a debugging dump of the matrix.

Parameters
----------
name : str
    Name of the dump."
