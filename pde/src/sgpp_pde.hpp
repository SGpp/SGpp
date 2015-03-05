// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PDE_HPP
#define PDE_HPP


#include "sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp"
#include "sgpp/pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichlet.hpp"
#include "sgpp/pde/application/HeatEquationSolver.hpp"
#include "sgpp/pde/application/HeatEquationSolverWithStretching.hpp"
#include "sgpp/pde/application/PoissonEquationSolver.hpp"
#include "sgpp/pde/operation/hash/OperationLaplaceLinear.hpp"
#include "sgpp/pde/operation/hash/OperationLaplaceLinearBoundary.hpp"
#include "sgpp/pde/operation/hash/OperationLaplaceLinearStretched.hpp"
#include "sgpp/pde/operation/hash/OperationLaplaceLinearStretchedBoundary.hpp"
#include "sgpp/pde/operation/hash/OperationLaplaceModLinear.hpp"
#include "sgpp/pde/operation/hash/OperationParabolicPDESolverSystemFreeBoundaries.hpp"
#include "sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitPeriodic.hpp"
#include "sgpp/pde/operation/hash/OperationMatrixLTwoDotPeriodic.hpp"

#include "sgpp/pde/operation/PdeOpFactory.hpp"

#endif /* PDE_HPP */