// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// Shared pointers.
%shared_ptr(sgpp::finance::BlackScholesParabolicPDESolverSystem)
%shared_ptr(sgpp::finance::BlackScholesParabolicPDESolverSystemEuroAmer)
%shared_ptr(sgpp::finance::BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP)

// The Good, i.e. without any modifications
%include "finance/src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystem.hpp"
%include "finance/src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp"
%include "finance/src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP.hpp"

%include "finance/src/sgpp/finance/application/BlackScholesSolver.hpp"
%include "finance/src/sgpp/finance/application/BlackScholesSolverWithStretching.hpp"

%include "finance/src/sgpp/finance/tools/VariableDiscountFactor.hpp"

%include "OpFactory.i"

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point }; 
