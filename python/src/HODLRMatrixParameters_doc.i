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
