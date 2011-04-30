/******************************************************************************
* Copyright (C) 201 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef SGPP_MPI_HPP
#define SGPP_MPI_HPP

#include "tools/MPI/SGppMPITools.hpp"

#include "solver/sle/ConjugateGradientsMPI.hpp"

#include "application/pde/PoissonEquationSolverMPI.hpp"
#include "application/pde/HeatEquationSolverMPI.hpp"

#include "algorithm/pde/PoissonEquationEllipticPDESolverSystemDirichletParallelMPI.hpp"
#include "algorithm/pde/HeatEquationParabolicPDESolverSystemParallelMPI.hpp"

#include "sgpp.hpp"

#endif /* SGPP_MPI_HPP */
