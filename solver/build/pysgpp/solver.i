// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%rename(SolverModuleSLESolver)          sgpp::solver::SLESolver;
%rename(SolverModuleBiCGStab)           sgpp::solver::BiCGStab;

// The Good, i.e. without any modifications
%include "solver/src/sgpp/solver/SGSolver.hpp"
%include "solver/src/sgpp/solver/SLESolver.hpp"
%include "solver/src/sgpp/solver/ODESolver.hpp"
%feature("director") ConjugateGradients;
%include "solver/src/sgpp/solver/sle/ConjugateGradients.hpp"
%include "solver/src/sgpp/solver/sle/BiCGStab.hpp"
%include "solver/src/sgpp/solver/ode/Euler.hpp"
%include "solver/src/sgpp/solver/ode/CrankNicolson.hpp"
%include "solver/src/sgpp/solver/TypesSolver.hpp"
%include "solver/src/sgpp/solver/SLESolverTypeParser.hpp"

%include "solver/src/sgpp/solver/operation/hash/OperationParabolicPDESolverSystem.hpp"

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point }; 
