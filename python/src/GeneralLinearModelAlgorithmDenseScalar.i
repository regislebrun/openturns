// SWIG file GeneralLinearModelAlgorithmDenseScalar.i

%{
#include "openturns/GeneralLinearModelAlgorithmDenseScalar.hxx"
%}

%include GeneralLinearModelAlgorithmDenseScalar_doc.i

%include openturns/GeneralLinearModelAlgorithmDenseScalar.hxx

namespace OT{ %extend GeneralLinearModelAlgorithmDenseScalar { GeneralLinearModelAlgorithmDenseScalar(const GeneralLinearModelAlgorithmDenseScalar & other) { return new OT::GeneralLinearModelAlgorithmDenseScalar(other); } } }

