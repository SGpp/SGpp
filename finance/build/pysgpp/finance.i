/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Joerg Blank (blankj@in.tum.de), Alexander Heinecke (alexander.heinecke@mytum.de)

// The Good, i.e. without any modifications
%include "finance/src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystem.hpp"
%include "finance/src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp"
%include "finance/src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP.hpp"

%include "finance/src/sgpp/finance/application/BlackScholesSolver.hpp"
%include "finance/src/sgpp/finance/application/BlackScholesSolverWithStretching.hpp"

%include "finance/src/sgpp/finance/tools/VariableDiscountFactor.hpp"

%include "finance/src/sgpp/finance/operation/FinanceOpFactory.hpp"

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%apply std::vector<std::pair<size_t, float_t> > *OUTPUT { std::vector<std::pair<size_t, float_t> >& result };
%apply std::vector<float_t> *INPUT { std::vector<float_t>& point }; 
