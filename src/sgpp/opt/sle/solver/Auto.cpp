/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/sle/solver/Auto.hpp"
#include "opt/sle/solver/Armadillo.hpp"
#include "opt/sle/solver/BiCGStab.hpp"
#include "opt/sle/solver/Eigen.hpp"
#include "opt/sle/solver/Gmmpp.hpp"
#include "opt/sle/solver/UMFPACK.hpp"
#include "opt/tools/Printer.hpp"

#include <cstddef>
#include <algorithm>
#include <map>

namespace sg
{
namespace opt
{
namespace sle
{
namespace solver
{

const double Auto::MAX_NNZ_RATIO_FOR_SPARSE = 0.1;
const double Auto::MAX_NNZ_RATIO_FOR_GMMPP = 0.05;
const double Auto::ESTIMATE_NNZ_ROWS_SAMPLE_SIZE = 0.05;

/**
 * Add a solver to the vector of solvers, if the solver is supported
 * (e.g. SG++ configured and compiled for use with the solve) and it's not already
 * in the vector.
 * 
 * @param solver    linear solver
 * @param solvers   vector of solvers
 * @param supports  map indicating which solvers are supported
 */
void addSolver(Solver *solver, std::vector<Solver *> &solvers,
               const std::map<Solver *, bool> &supports)
{
    // add solver if it's supported and not already in the vector
    if ((supports.at(solver)) &&
        (std::find(solvers.begin(), solvers.end(), solver) == solvers.end()))
    {
        solvers.push_back(solver);
    }
}

bool Auto::solve(system::System &system, const std::vector<double> &b,
                 std::vector<double> &x) const
{
    std::vector<std::vector<double> > B;
    std::vector<std::vector<double> > X;
    B.push_back(b);
    X.push_back(x);
    
    if (solve(system, B, X))
    {
        x = X[0];
        return true;
    } else
    {
        return false;
    }
}

bool Auto::solve(system::System &system, const std::vector<std::vector<double> > &B,
                 std::vector<std::vector<double> > &X) const
{
    tools::printer.printStatusBegin("Solving linear system (automatic method)...");
    
    Armadillo solver_armadillo;
    Eigen solver_eigen;
    UMFPACK solver_umfpack;
    Gmmpp solver_gmmpp;
    BiCGStab solver_bicgstab;
    
    std::map<Solver *, bool> supports;
    
    // by default, only BiCGStab is supported
    supports[&solver_armadillo] = false;
    supports[&solver_eigen] = false;
    supports[&solver_umfpack] = false;
    supports[&solver_gmmpp] = false;
    supports[&solver_bicgstab] = true;
    
#ifdef USEARMADILLO
    supports[&solver_armadillo] = true;
#endif
    
#ifdef USEEIGEN
    supports[&solver_eigen] = true;
#endif
    
#ifdef USEUMFPACK
    supports[&solver_umfpack] = true;
#endif
    
#ifdef USEGMMPP
    supports[&solver_gmmpp] = true;
#endif
    
    // solvers to be used, the solver which should be tried first should be the first element
    std::vector<Solver *> solvers;
    const size_t n = system.getDimension();
    
    if (supports[&solver_umfpack] || supports[&solver_gmmpp])
    {
        // if at least one of the sparse solvers is supported
        // ==> estimate sparsity ratio of matrix by considering every inc-th row
        size_t nrows = 0;
        size_t nnz = 0;
        size_t inc = static_cast<size_t>(
                ESTIMATE_NNZ_ROWS_SAMPLE_SIZE * static_cast<double>(n)) + 1;
        
        tools::printer.printStatusUpdate("estimating sparsity pattern");
        
        for (size_t i = 0; i < n; i += inc)
        {
            nrows++;
            
            for (size_t j = 0; j < n; j++)
            {
                if (system.isMatrixEntryNonZero(i, j))
                {
                    nnz++;
                }
            }
        }
        
        // calculate estimate ratio nonzero entries
        double nnz_ratio = static_cast<double>(nnz) /
                (static_cast<double>(nrows) * static_cast<double>(n));
        
        // print ratio
        {
            char str[10];
            snprintf(str, 10, "%.1f%%", nnz_ratio * 100.0);
            tools::printer.printStatusUpdate("estimated nnz ratio: " + std::string(str));
            tools::printer.printStatusNewLine();
        }
        
        if ((n > MAX_DIM_FOR_FULL) || (nnz_ratio <= MAX_NNZ_RATIO_FOR_GMMPP))
        {
            if (B.size() == 1)
            {
                // prefer Gmm++ over UMFPACK
                addSolver(&solver_gmmpp, solvers, supports);
                addSolver(&solver_umfpack, solvers, supports);
            } else
            {
                // prefer UMFPACK over Gmm++ (UMFPACK can solve multiple systems simultaneously)
                addSolver(&solver_umfpack, solvers, supports);
                addSolver(&solver_gmmpp, solvers, supports);
            }
        } else if (nnz_ratio <= MAX_NNZ_RATIO_FOR_SPARSE)
        {
            // prefer UMFPACK over Gmm++
            addSolver(&solver_umfpack, solvers, supports);
            addSolver(&solver_gmmpp, solvers, supports);
        }
    }
    
    // add all remaining solvers (prefer Armadillo over Eigen)
    addSolver(&solver_armadillo, solvers, supports);
    addSolver(&solver_eigen, solvers, supports);
    addSolver(&solver_gmmpp, solvers, supports);
    addSolver(&solver_umfpack, solvers, supports);
    addSolver(&solver_bicgstab, solvers, supports);
    
    for (size_t i = 0; i < solvers.size(); i++)
    {
        // try solver
        bool result = solvers[i]->solve(system, B, X);
        
        if (result)
        {
            tools::printer.printStatusEnd();
            return true;
        } else if ((solvers[i] == &solver_gmmpp) && (n > MAX_DIM_FOR_FULL))
        {
            // don't use full solvers and return approximative solution
            tools::printer.printStatusEnd("warning: using non-converged solution of iterative "
                                          "solver, residual can be large "
                                          "(matrix too large to try other solvers)");
            return true;
        }
    }
    
    tools::printer.printStatusEnd("error: could not solve linear system!");
    return false;
}

}
}
}
}
