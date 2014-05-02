#include "opt/sle/SolverAuto.hpp"
#include "opt/sle/SolverArmadillo.hpp"
#include "opt/sle/SolverBiCGStab.hpp"
#include "opt/sle/SolverEigen.hpp"
#include "opt/sle/SolverUMFPACK.hpp"
#include "opt/tools/Printer.hpp"

#include <cstddef>
#include <algorithm>

namespace sg
{
namespace opt
{
namespace sle
{

const double SolverAuto::MAX_NNZ_RATIO_FOR_SPARSE = 0.2;
const double SolverAuto::ESTIMATE_NNZ_ROWS_SAMPLE_SIZE = 0.05;

bool SolverAuto::solve(System &system, std::vector<double> &x) const
{
    tools::printer.printStatusBegin("Solving linear system (automatic method)...");
    
    bool armadillo_support = false;
    bool eigen_support = false;
    bool umfpack_support = false;
    
#ifdef USEARMADILLO
    armadillo_support = true;
#endif
    
#ifdef USEEIGEN
    eigen_support = true;
#endif
    
#ifdef USEUMFPACK
    umfpack_support = true;
#endif
    
    bool solve_with_armadillo = false;
    bool solve_with_eigen = false;
    bool solve_with_umfpack = false;
    
    size_t n = system.getDimension();
    
    /*const std::vector<double> &b = system.getRHS();
    size_t nnz = 0;
    std::vector<uint32_t> Ti;
    std::vector<uint32_t> Tj;
    std::vector<double> Tx;*/
    //bool sparse_already_constructed = false;
    
    if (umfpack_support)
    {
        /*for (uint32_t i = 0; i < n; i++)
        {
            if (i % 100 == 0)
            {
                char str[10];
                snprintf(str, 10, "%.1f%%",
                         static_cast<double>(i) / static_cast<double>(n) * 100.0);
                tools::printer.printStatusUpdate("constructing sparse matrix (" +
                                                 std::string(str) + ")");
            }
            
            for (uint32_t j = 0; j < n; j++)
            {
                double entry = system.getMatrixEntry(i, j);
                
                if (entry != 0.0)
                {
                    Ti.push_back(i);
                    Tj.push_back(j);
                    Tx.push_back(entry);
                    nnz++;
                }
            }
        }
        
        tools::printer.printStatusNewLine();
        tools::printer.printStatusUpdate("nnz = " + std::to_string(nnz));
        tools::printer.printStatusNewLine();
        sparse_already_constructed = true;
        
        if (static_cast<double>(nnz) <= MAX_NNZ_RATIO_FOR_SPARSE *
                                        static_cast<double>(n)*static_cast<double>(n))
        {
            solve_with_umfpack = true;
        }*/
        
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
        
        /*tools::printer.printStatusUpdate("estimating sparsity pattern");
        
        size_t nnz = 0;
        size_t nrows = 1 + static_cast<int>(ESTIMATE_NNZ_ROWS_SAMPLE_SIZE *
                                            static_cast<double>(n));
        std::vector<size_t> row_ind(n, 0.0);
        
        for (size_t k = 0; k < n; k++)
        {
            row_ind[k] = k;
        }
        
        // caution: shuffle is not random, since srand is not called at this time
        std::random_shuffle(row_ind.begin(), row_ind.end());
        
        for (size_t k = 0; k < nrows; k++)
        {
            size_t i = row_ind[k];
            
            for (size_t j = 0; j < n; j++)
            {
                if (system.isMatrixEntryNonZero(i, j))
                {
                    nnz++;
                }
            }
        }*/
        
        double nnz_ratio = static_cast<double>(nnz) /
                (static_cast<double>(nrows) * static_cast<double>(n));
        
        {
            char str[10];
            snprintf(str, 10, "%.1f%%", nnz_ratio * 100.0);
            tools::printer.printStatusUpdate("estimated nnz ratio: " + std::string(str));
            tools::printer.printStatusNewLine();
        }
        
        if (nnz_ratio <= MAX_NNZ_RATIO_FOR_SPARSE)
        {
            solve_with_umfpack = true;
        }
    }
    
    if (!solve_with_umfpack)
    {
        if (armadillo_support)
        {
            solve_with_armadillo = true;
        } else if (eigen_support)
        {
            solve_with_eigen = true;
        }
    }
    
    bool result;
    
    if (solve_with_armadillo)
    {
        SolverArmadillo solver;
        
        /*if (sparse_already_constructed)
        {
            result = solver.solve(Ti, Tj, Tx, b, x);
        } else
        {*/
            result = solver.solve(system, x);
        //}
    } else if (solve_with_eigen)
    {
        SolverEigen solver;
        
        /*if (sparse_already_constructed)
        {
            result = solver.solve(Ti, Tj, Tx, b, x);
        } else
        {*/
            result = solver.solve(system, x);
        //}
    } else if (solve_with_umfpack)
    {
        SolverUMFPACK solver;
        
        /*if (sparse_already_constructed)
        {
            result = solver.solve(Ti, Tj, Tx, b, x);
        } else
        {*/
            result = solver.solve(system, x);
        //}
    } else
    {
        SolverBiCGStab solver;
        
        /*if (sparse_already_constructed)
        {
            result = solver.solve(Ti, Tj, Tx, b, x);
        } else
        {*/
            result = solver.solve(system, x);
        //}
    }
    
    tools::printer.printStatusEnd();
    
    return result;
}

/*bool SolverAuto::solve(const std::vector<uint32_t> &Ti, const std::vector<uint32_t> &Tj,
                       const std::vector<double> &Tx, const std::vector<double> &b,
                       std::vector<double> &x) const
{
    tools::printer.printStatusBegin("Solving sparse linear system (automatic method)...");
    
    bool armadillo_support = false;
    bool eigen_support = false;
    bool umfpack_support = false;
    
#ifdef USEARMADILLO
    armadillo_support = true;
#endif
    
#ifdef USEEIGEN
    eigen_support = true;
#endif
    
#ifdef USEUMFPACK
    umfpack_support = true;
#endif
    
    bool solve_with_armadillo = false;
    bool solve_with_eigen = false;
    bool solve_with_umfpack = false;
    
    if (umfpack_support)
    {
        double nnz = static_cast<double>(Tx.size());
        double n = static_cast<double>(b.size());
        
        if (static_cast<double>(nnz) <= MAX_NNZ_RATIO_FOR_SPARSE *
                                        static_cast<double>(n)*static_cast<double>(n))
        {
            solve_with_umfpack = true;
        }
    }
    
    if (!solve_with_umfpack)
    {
        if (armadillo_support)
        {
            solve_with_armadillo = true;
        } else if (eigen_support)
        {
            solve_with_eigen = true;
        }
    }
    
    bool result;
    
    if (solve_with_armadillo)
    {
        SolverArmadillo solver;
        result = solver.solve(Ti, Tj, Tx, b, x);
    } else if (solve_with_eigen)
    {
        SolverEigen solver;
        result = solver.solve(Ti, Tj, Tx, b, x);
    } else if (solve_with_umfpack)
    {
        SolverUMFPACK solver;
        result = solver.solve(Ti, Tj, Tx, b, x);
    } else
    {
        SolverBiCGStab solver;
        result = solver.solve(Ti, Tj, Tx, b, x);
    }
    
    tools::printer.printStatusEnd();
    
    return result;
}*/

}
}
}
