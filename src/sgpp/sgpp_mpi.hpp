/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef SGPP_MPI_HPP
#define SGPP_MPI_HPP

#include "parallel/parallel/tools/MPI/SGppMPITools.hpp"

#include "solver/sle/ConjugateGradientsMPI.hpp"
#include "solver/sle/BiCGStabMPI.hpp"

#include "pde/application/PoissonEquationSolverMPI.hpp"
#include "pde/application/HeatEquationSolverMPI.hpp"
#include "pde/application/BlackScholesSolverMPI.hpp"

#include "pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichletParallelMPI.hpp"
#include "pde/algorithm/HeatEquationParabolicPDESolverSystemParallelMPI.hpp"
#include "pde/algorithm/BlackScholesParabolicPDESolverSystemEuroAmerParallelMPI.hpp"

#endif /* SGPP_MPI_HPP */
