// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// Shared pointers.
%shared_ptr(sgpp::pde::UpDownOneOpDim)
%shared_ptr(sgpp::pde::UpDownOneOpDimWithShadow)
%shared_ptr(sgpp::pde::OperationMatrixLTwoDotExplicitPeriodic)
%shared_ptr(sgpp::pde::OperationMatrixLTwoDotPeriodic)
%shared_ptr(sgpp::pde::OperationParabolicPDESolverSystemDirichlet)
%shared_ptr(sgpp::pde::HeatEquationParabolicPDESolverSystem)
%shared_ptr(sgpp::pde::OperationParabolicPDESolverSystemFreeBoundaries)
%shared_ptr(sgpp::solver::OperationParabolicPDESolverSystem)
%shared_ptr(sgpp::pde::OperationLaplaceLinear)
%shared_ptr(sgpp::pde::OperationLaplaceLinearBoundary)
%shared_ptr(sgpp::pde::OperationLaplaceModLinear)
%shared_ptr(sgpp::pde::OperationLaplacePrewavelet)
%shared_ptr(sgpp::pde::OperationLaplaceLinearStretched)
%shared_ptr(sgpp::pde::OperationLaplaceLinearStretchedBoundary)

// The Good, i.e. without any modifications
%include "solver/src/sgpp/solver/operation/hash/OperationParabolicPDESolverSystem.hpp"
%include "pde/src/sgpp/pde/operation/hash/OperationParabolicPDESolverSystemDirichlet.hpp"

%include "pde/src/sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp"

%include "pde/src/sgpp/pde/application/PDESolver.hpp"
%include "pde/src/sgpp/pde/application/ParabolicPDESolver.hpp"
%include "pde/src/sgpp/pde/application/HeatEquationSolver.hpp"
%include "pde/src/sgpp/pde/application/HeatEquationSolver.hpp"
%include "pde/src/sgpp/pde/application/HeatEquationSolverWithStretching.hpp"
%include "pde/src/sgpp/pde/application/EllipticPDESolver.hpp"
%include "pde/src/sgpp/pde/application/PoissonEquationSolver.hpp"

%include "OpFactory.i"

%include "pde/src/sgpp/pde/operation/hash/OperationParabolicPDESolverSystemFreeBoundaries.hpp"
%include "pde/src/sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitPeriodic.hpp"
%include "pde/src/sgpp/pde/operation/hash/OperationMatrixLTwoDotPeriodic.hpp"

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point };
