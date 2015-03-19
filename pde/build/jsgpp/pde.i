// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// The Good, i.e. without any modifications
%include "pde/src/sgpp/pde/operation/hash/OperationParabolicPDESolverSystem.hpp"
%include "pde/src/sgpp/pde/operation/hash/OperationParabolicPDESolverSystemDirichlet.hpp"

%include "pde/src/sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp"

%include "pde/src/sgpp/pde/application/PDESolver.hpp"
%include "pde/src/sgpp/pde/application/ParabolicPDESolver.hpp"
%include "pde/src/sgpp/pde/application/HeatEquationSolver.hpp"
%include "pde/src/sgpp/pde/application/HeatEquationSolver.hpp"
%include "pde/src/sgpp/pde/application/HeatEquationSolverWithStretching.hpp"
%include "pde/src/sgpp/pde/application/EllipticPDESolver.hpp"
%include "pde/src/sgpp/pde/application/PoissonEquationSolver.hpp"

%include "pde/src/sgpp/pde/operation/hash/OperationParabolicPDESolverSystemFreeBoundaries.hpp"
%include "pde/src/sgpp/pde/operation/PdeOpFactory.hpp"
%include "pde/src/sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitPeriodic.hpp"
%include "pde/src/sgpp/pde/operation/hash/OperationMatrixLTwoDotPeriodic.hpp"

//%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

//%apply std::vector<std::pair<size_t, float_t> > *OUTPUT { std::vector<std::pair<size_t, float_t> >& result };
//%apply std::vector<float_t> *INPUT { std::vector<float_t>& point }; 
