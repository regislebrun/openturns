%feature("docstring") OT::HODLRMatrixParameters
"Parameters for HODLRMatrix class.

This class regroups the parameters used by :class:`~openturns.HODLRMatrix`."

%feature("docstring") OT::HODLRMatrixParameters::getAssemblyEpsilon
"Return the assembly epsilon.

Returns
-------
epsilon : float
    Assembly epsilon for ACA compression."

%feature("docstring") OT::HODLRMatrixParameters::setAssemblyEpsilon
"Set the assembly epsilon.

Parameters
----------
epsilon : float
    Assembly epsilon for ACA compression."

%feature("docstring") OT::HODLRMatrixParameters::getRecompressionEpsilon
"Return the recompression epsilon.

Returns
-------
epsilon : float
    Recompression epsilon."

%feature("docstring") OT::HODLRMatrixParameters::setRecompressionEpsilon
"Set the recompression epsilon.

Parameters
----------
epsilon : float
    Recompression epsilon."

%feature("docstring") OT::HODLRMatrixParameters::getMinLeafSize
"Return the minimum leaf size.

Returns
-------
size : int
    Minimum number of rows/columns for a leaf block."

%feature("docstring") OT::HODLRMatrixParameters::setMinLeafSize
"Set the minimum leaf size.

Parameters
----------
size : int
    Minimum number of rows/columns for a leaf block."

%feature("docstring") OT::HODLRMatrixParameters::getMaxRank
"Return the maximum rank for low-rank blocks.

Returns
-------
rank : int
    Maximum rank of the low-rank blocks. Zero means the rank is
    adaptive, i.e. driven by the assembly epsilon."

%feature("docstring") OT::HODLRMatrixParameters::setMaxRank
"Set the maximum rank for low-rank blocks.

Parameters
----------
rank : int
    Maximum rank of the low-rank blocks. Zero (the default) means the
    rank is adaptive, i.e. each block is compressed up to the assembly
    epsilon. A positive value caps the rank of every block; blocks that
    hit the cap before reaching the assembly epsilon are reported by a
    warning as rank-starved."

%feature("docstring") OT::HODLRMatrixParameters::getFactorizationMethod
"Return the factorization method.

Returns
-------
method : str
    Factorization method, either 'LU' or 'LLt'."

%feature("docstring") OT::HODLRMatrixParameters::setFactorizationMethod
"Set the factorization method.

Parameters
----------
method : str
    Factorization method, either 'LU' or 'LLt'."

%feature("docstring") OT::HODLRMatrixParameters::getMaxRank
"Return the maximum rank of the low-rank blocks.

Returns
-------
maxRank : int
    Maximum rank used for the low-rank approximation during assembly."

%feature("docstring") OT::HODLRMatrixParameters::setMaxRank
"Set the maximum rank of the low-rank blocks.

Parameters
----------
maxRank : int
    Maximum rank used for the low-rank approximation during assembly."

%feature("docstring") OT::HODLRMatrixParameters::getUseSpatialOrdering
"Return whether the spatial ordering is used.

Returns
-------
use : bool
    True if the vertices are reordered along a space-filling curve
    before assembly."

%feature("docstring") OT::HODLRMatrixParameters::setUseSpatialOrdering
"Set whether the spatial ordering is used.

Parameters
----------
use : bool
    If True (the default), the vertices are reordered along a
    space-filling curve before assembly, so that the recursive split
    of the HODLR tree separates spatially close points at the leaves."
