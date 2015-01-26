// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef SGPP_MPI_HPP
#define SGPP_MPI_HPP

#include "parallel/tools/MPI/SGppMPITools.hpp"

#include "parallel/solver/sle/ConjugateGradientsMPI.hpp"
#include "parallel/solver/sle/BiCGStabMPI.hpp"

#include "parallel/pde/application/PoissonEquationSolverMPI.hpp"
#include "parallel/pde/application/HeatEquationSolverMPI.hpp"
#include "parallel/finance/application/BlackScholesSolverMPI.hpp"

#include "parallel/pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichletParallelMPI.hpp"
#include "parallel/pde/algorithm/HeatEquationParabolicPDESolverSystemParallelMPI.hpp"
#include "parallel/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmerParallelMPI.hpp"

#endif /* SGPP_MPI_HPP */