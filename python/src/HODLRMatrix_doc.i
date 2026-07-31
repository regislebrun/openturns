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
    Factorization method, either 'LU' or 'LLt'."

%feature("docstring") OT::HODLRMatrix::solve
"Solve a linear system.

Parameters
----------
b : 1-d or 2-d sequence of float
    Right-hand side.
trans : bool, optional
    Whether to solve with the transpose. Default is False.

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
