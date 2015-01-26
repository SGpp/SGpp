/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de

#ifndef PDE_HPP
#define PDE_HPP


#include "pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp"
#include "pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichlet.hpp"
#include "pde/application/HeatEquationSolver.hpp"
#include "pde/application/HeatEquationSolverWithStretching.hpp"
#include "pde/application/PoissonEquationSolver.hpp"
#include "pde/basis/linear/noboundary/operation/OperationLaplaceLinear.hpp"
#include "pde/basis/linear/boundary/operation/OperationLaplaceLinearBoundary.hpp"
#include "pde/basis/linearstretched/noboundary/operation/OperationLaplaceLinearStretched.hpp"
#include "pde/basis/linearstretched/boundary/operation/OperationLaplaceLinearStretchedBoundary.hpp"
#include "pde/basis/modlinear/operation/OperationLaplaceModLinear.hpp"
#include "pde/operation/OperationParabolicPDESolverSystemFreeBoundaries.hpp""
#include "pde/basis/periodic/operation/OperationMatrixLTwoDotExplicitPeriodic.hpp"
#include "pde/basis/periodic/operation/OperationMatrixLTwoDotPeriodic.hpp"

#include "pde/operation/PdeOpFactory.hpp"

#endif /* PDE_HPP */
