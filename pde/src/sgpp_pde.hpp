/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de

#ifndef PDE_HPP
#define PDE_HPP


#include "sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp"
#include "sgpp/pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichlet.hpp"
#include "sgpp/pde/application/HeatEquationSolver.hpp"
#include "sgpp/pde/application/HeatEquationSolverWithStretching.hpp"
#include "sgpp/pde/application/PoissonEquationSolver.hpp"
#include "sgpp/pde/basis/linear/noboundary/operation/OperationLaplaceLinear.hpp"
#include "sgpp/pde/basis/linear/boundary/operation/OperationLaplaceLinearBoundary.hpp"
#include "sgpp/pde/basis/linearstretched/noboundary/operation/OperationLaplaceLinearStretched.hpp"
#include "sgpp/pde/basis/linearstretched/boundary/operation/OperationLaplaceLinearStretchedBoundary.hpp"
#include "sgpp/pde/basis/modlinear/operation/OperationLaplaceModLinear.hpp"
#include "sgpp/pde/operation/OperationParabolicPDESolverSystemFreeBoundaries.hpp"
#include "sgpp/pde/basis/periodic/operation/OperationMatrixLTwoDotExplicitPeriodic.hpp"
#include "sgpp/pde/basis/periodic/operation/OperationMatrixLTwoDotPeriodic.hpp"

#include "sgpp/pde/operation/PdeOpFactory.hpp"

#endif /* PDE_HPP */
