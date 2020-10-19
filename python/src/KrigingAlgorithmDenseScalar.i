// SWIG file KrigingAlgorithmDenseScalar.i

%{
#include "openturns/KrigingAlgorithmDenseScalar.hxx"
%}

%include KrigingAlgorithmDenseScalar_doc.i

%include openturns/KrigingAlgorithmDenseScalar.hxx

namespace OT{ %extend KrigingAlgorithmDenseScalar { KrigingAlgorithmDenseScalar(const KrigingAlgorithmDenseScalar & other) { return new OT::KrigingAlgorithmDenseScalar(other); } } }
