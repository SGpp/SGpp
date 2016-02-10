// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// %include "combigrid/src/sgpp/combigrid/utils/combigrid_ultils.hpp"
// %ignore combigrid::CombigridLevelVector::operator=;
// %include "combigrid/src/sgpp/combigrid/utils/CombigridLevelVector.hpp"  
// %include "combigrid/src/sgpp/combigrid/basisfunction/CombiBasisFunctionBasis.hpp"
// %include "combigrid/src/sgpp/combigrid/basisfunction/CombiLinearBasisFunction.hpp"
// %include "combigrid/src/sgpp/combigrid/domain/AbstractStretchingMaker.hpp"
// %include "combigrid/src/sgpp/combigrid/domain/CombiDomain1D.hpp" 
// %include "combigrid/src/sgpp/combigrid/domain/CombiGridDomain.hpp"
// %include "combigrid/src/sgpp/combigrid/domain/CombiAtanSpecialStretching.hpp"
// %include "combigrid/src/sgpp/combigrid/domain/CombiTanStretching.hpp"
// %include "combigrid/src/sgpp/combigrid/combischeme/CombiTS_CT.hpp"
// %include "combigrid/src/sgpp/combigrid/combischeme/CombiS_CT.hpp"
// %include "combigrid/src/sgpp/combigrid/combischeme/CombiArbitraryScheme.hpp"
// %include "combigrid/src/sgpp/combigrid/combigrid/SerialCombiGrid.hpp"
// 
// %rename(__add__) combigrid::CombigridLevelVector::operator+;
// %rename(__mul__) combigrid::CombigridLevelVector::operator*;
// %rename(__sub__) combigrid::CombigridLevelVector::operator-;
// %rename(__new__) combigrid::CombigridLevelVector::operator=; 
// 
// %include "combigrid/src/sgpp/combigrid/fullgrid/CombiFullGrid.hpp"
// %template(doubleFullGrid) combigrid::FullGrid<double>;
// 
%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%apply std::vector<std::pair<size_t, SGPP::float_t> > *OUTPUT { std::vector<std::pair<size_t, SGPP::float_t> >& result };
%apply std::vector<SGPP::float_t> *INPUT { std::vector<SGPP::float_t>& point }; 
