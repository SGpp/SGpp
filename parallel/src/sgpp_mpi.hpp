// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_MPI_HPP
#define SGPP_MPI_HPP

#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>

#include <sgpp/parallel/solver/sle/ConjugateGradientsMPI.hpp>
#include <sgpp/parallel/solver/sle/BiCGStabMPI.hpp>

#include <sgpp/parallel/pde/application/PoissonEquationSolverMPI.hpp>
#include <sgpp/parallel/pde/application/HeatEquationSolverMPI.hpp>
#include <sgpp/parallel/finance/application/BlackScholesSolverMPI.hpp>

#include <sgpp/parallel/pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichletParallelMPI.hpp>
#include <sgpp/parallel/pde/algorithm/HeatEquationParabolicPDESolverSystemParallelMPI.hpp>
#include <sgpp/parallel/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmerParallelMPI.hpp>

#endif /* SGPP_MPI_HPP */
