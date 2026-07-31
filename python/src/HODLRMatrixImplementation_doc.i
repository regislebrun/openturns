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
