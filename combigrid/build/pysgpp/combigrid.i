/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Joerg Blank (blankj@in.tum.de), Alexander Heinecke (alexander.heinecke@mytum.de)

%module(directors="1") combigrid
%feature("docstring");

%include "../../../base/src/sgpp/globaldef.hpp"

%include "stl.i"
%include "std_vector.i"
%include "std_pair.i"
%include "std_complex.i"
%include "std_map.i"

%include "cpointer.i" 
%include "typemaps.i"

%include "exception.i"

%{
#define SWIG_FILE_WITH_INIT
%}

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%include "carrays.i"
%array_class(unsigned int, unsignedIntArray);
%array_class(bool,BoolArray);
%array_class(int, IntArray);


namespace std {
	%template(IntVector) vector<int>;
	%template(IntIntVector) vector< vector<int> >; 
	%template(BoolVector) vector<bool>;
	%template(DoubleVector) vector<double>;
	%template(IndexValPair) pair<size_t, double>;
        %template(IndexValVector) vector<pair<size_t, double> >;
        // For OnlinePredictiveRefinementDimension
        %template(refinement_key) std::pair<size_t, unsigned int>;
        %template(refinement_map) std::map<std::pair<size_t, unsigned int>, double>;

}

// This should include all necessary header files
%{
#include "src/combigrid.hpp"
%}

// The Good, i.e. without any modifications
//%include "FullGrid.i"
//%include "combigrid/src/sgpp/base/grid/combination/FullGridSet.hpp"
//%include "FullGridSet.i"

%include "combigrid/src/sgpp/combigrid/utils/combigrid_ultils.hpp"
%ignore combigrid::CombigridLevelVector::operator=;
%include "combigrid/src/sgpp/combigrid/utils/CombigridLevelVector.hpp"  
%include "combigrid/src/sgpp/combigrid/basisfunction/CombiBasisFunctionBasis.hpp"
%include "combigrid/src/sgpp/combigrid/basisfunction/CombiLinearBasisFunction.hpp"
%include "combigrid/src/sgpp/combigrid/domain/AbstractStretchingMaker.hpp"
%include "combigrid/src/sgpp/combigrid/domain/CombiDomain1D.hpp" 
%include "combigrid/src/sgpp/combigrid/domain/CombiGridDomain.hpp"
%include "combigrid/src/sgpp/combigrid/domain/CombiAtanSpecialStretching.hpp"
%include "combigrid/src/sgpp/combigrid/domain/CombiTanStretching.hpp"
%include "combigrid/src/sgpp/combigrid/domain/CombiUniformStretching.hpp"   
%include "combigrid/src/sgpp/combigrid/combischeme/CombiSchemeBasis.hpp" 
%include "combigrid/src/sgpp/combigrid/combischeme/CombiTS_CT.hpp"
%include "combigrid/src/sgpp/combigrid/combischeme/CombiS_CT.hpp"
%include "combigrid/src/sgpp/combigrid/combischeme/CombiArbitraryScheme.hpp"
%include "combigrid/src/sgpp/combigrid/combigridkernel/CombiGridKernel.hpp"
%include "combigrid/src/sgpp/combigrid/combigrid/AbstractCombiGrid.hpp"
%include "combigrid/src/sgpp/combigrid/combigrid/SerialCombiGrid.hpp"
%include "combigrid/src/sgpp/combigrid/combigrid/AdaptiveSerialCombiGrid.hpp"
%include "combigrid/src/sgpp/combigrid/combigrid/AdaptiveSerialCombiGridVariableCoefficients.hpp" 

%rename(__add__) combigrid::CombigridLevelVector::operator+;
%rename(__mul__) combigrid::CombigridLevelVector::operator*;
%rename(__sub__) combigrid::CombigridLevelVector::operator-;
%rename(__new__) combigrid::CombigridLevelVector::operator=; 

//%template(ComplexDouble) complex<double>;
//
%include "combigrid/src/sgpp/combigrid/fullgrid/CombiFullGrid.hpp"
%template(doubleFullGrid) combigrid::FullGrid<double>;

//%template(FullGridC) combigrid::FullGrid< complex<double> >;
//%template(CombiGridKernelC) combigrid::CombiGridKernel< complex<double> >;
%template(CombiGridKernelD) combigrid::CombiGridKernel< double >;   
//%template(ComplexVector) std::vector< complex<double> >;

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point }; 
