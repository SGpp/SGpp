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
const double Auto::ESTIMATE_NNZ_ROWS_SAMPLE_SIZE = 0.05;
const double Auto::MAX_NNZ_RATIO_FOR_GMMPP = 0.05;

void addSolver(Solver *solver, std::vector<Solver *> &solvers,
               const std::map<Solver *, bool> &supports)
{
    if ((supports.at(solver)) &&
        (std::find(solvers.begin(), solvers.end(), solver) == solvers.end()))
    {
        solvers.push_back(solver);
    }
}

bool Auto::solve(system::System &system, std::vector<double> &x) const
{
    tools::printer.printStatusBegin("Solving linear system (automatic method)...");
    
    Armadillo solver_armadillo;
    Eigen solver_eigen;
    UMFPACK solver_umfpack;
    Gmmpp solver_gmmpp;
    BiCGStab solver_bicgstab;
    
    std::map<Solver *, bool> supports =
    {
        {&solver_armadillo, false},
        {&solver_eigen, false},
        {&solver_umfpack, false},
        {&solver_gmmpp, false},
        {&solver_bicgstab, true}
    };
    
#ifdef USEARMADILLO
    //armadillo_support = true;
    supports[&solver_armadillo] = true;
#endif
    
#ifdef USEEIGEN
    //eigen_support = true;
    supports[&solver_eigen] = true;
#endif
    
#ifdef USEUMFPACK
    //umfpack_support = true;
    supports[&solver_umfpack] = true;
#endif
    
#ifdef USEGMMPP
    //gmmpp_support = true;
    supports[&solver_gmmpp] = true;
#endif
    
    std::vector<Solver *> solvers;
    size_t n = system.getDimension();
    
    if (supports[&solver_umfpack] || supports[&solver_gmmpp])
    {
        size_t nrows = 0;
        size_t nnz = 0;
        size_t inc = static_cast<int>(ESTIMATE_NNZ_ROWS_SAMPLE_SIZE * static_cast<double>(n)) + 1;
        
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
        
        double nnz_ratio = static_cast<double>(nnz) /
                (static_cast<double>(nrows) * static_cast<double>(n));
        
        {
            char str[10];
            snprintf(str, 10, "%.1f%%", nnz_ratio * 100.0);
            tools::printer.printStatusUpdate("estimated nnz ratio: " + std::string(str));
            tools::printer.printStatusNewLine();
        }
        
        if ((n > MAX_DIM_FOR_FULL) || (nnz_ratio <= MAX_NNZ_RATIO_FOR_GMMPP))
        {
            addSolver(&solver_gmmpp, solvers, supports);
            addSolver(&solver_umfpack, solvers, supports);
        } else if (nnz_ratio <= MAX_NNZ_RATIO_FOR_SPARSE)
        {
            addSolver(&solver_umfpack, solvers, supports);
            addSolver(&solver_gmmpp, solvers, supports);
        }
    }
    
    addSolver(&solver_armadillo, solvers, supports);
    addSolver(&solver_eigen, solvers, supports);
    addSolver(&solver_gmmpp, solvers, supports);
    addSolver(&solver_umfpack, solvers, supports);
    addSolver(&solver_bicgstab, solvers, supports);
    
    for (size_t i = 0; i < solvers.size(); i++)
    {
        bool result = solvers[i]->solve(system, x);
        
        if (result)
        {
            tools::printer.printStatusEnd();
            return true;
        } else if ((solvers[i] == &solver_gmmpp) && (n > MAX_DIM_FOR_FULL))
        {
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
