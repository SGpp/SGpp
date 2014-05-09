#include "opt/sle/solver/UMFPACK.hpp"
#include "opt/sle/system/Parallelizable.hpp"
#include "opt/tools/Printer.hpp"

#ifdef USEUMFPACK
#include <suitesparse/umfpack.h>
#endif

#include <cstddef>
#include <iostream>
#include <algorithm>

namespace sg
{
namespace opt
{
namespace sle
{
namespace solver
{

bool UMFPACK::solve(system::System &system, std::vector<double> &x) const
{
#ifdef USEUMFPACK
    tools::printer.printStatusBegin("Solving linear system (UMFPACK)...");
    
    size_t n = system.getDimension();
    const std::vector<double> &b = system.getRHS();
    
    size_t nnz = 0;
    std::vector<uint32_t> Ti;
    std::vector<uint32_t> Tj;
    std::vector<double> Tx;
    
    #pragma omp parallel if (system.isParallelizable())
    {
        system::System *real_system;
        
        if (system.isParallelizable())
        {
            real_system = dynamic_cast<system::Parallelizable &>(system).clone();
        } else
        {
            real_system = &system;
        }
        
        #pragma omp for ordered schedule(dynamic)
        for (uint32_t i = 0; i < n; i++)
        {
            for (uint32_t j = 0; j < n; j++)
            {
                double entry = real_system->getMatrixEntry(i, j);
                
                if (entry != 0)
                {
                    #pragma omp critical
                    {
                        Ti.push_back(i);
                        Tj.push_back(j);
                        Tx.push_back(entry);
                        nnz++;
                    }
                }
            }
            
            if (i % 100 == 0)
            {
                #pragma omp ordered
                {
                    char str[10];
                    snprintf(str, 10, "%.1f%%",
                             static_cast<double>(i) / static_cast<double>(n) * 100.0);
                    tools::printer.printStatusUpdate("constructing sparse matrix (" +
                                                     std::string(str) + ")");
                }
            }
        }
        
        if (system.isParallelizable())
        {
            delete real_system;
        }
    }
    
    // TODO: remove this if not needed anymore
    tools::printer.printVectorToFile("data/Ti.dat", Ti);
    tools::printer.printVectorToFile("data/Tj.dat", Tj);
    tools::printer.printVectorToFile("data/Tx.dat", Tx);
    
    tools::printer.printStatusUpdate("constructing sparse matrix (100.0%)");
    tools::printer.printStatusNewLine();
    
    {
        char str[10];
        double nnz_ratio = static_cast<double>(nnz) /
                           (static_cast<double>(n) * static_cast<double>(n));
        snprintf(str, 10, "%.1f%%", nnz_ratio * 100.0);
        tools::printer.printStatusUpdate("nnz ratio: " + std::string(str));
        tools::printer.printStatusNewLine();
    }
    
    bool result = solve(Ti, Tj, Tx, b, x, false);
    
    return result;
#else
    std::cerr << "Error in sg::opt::sle::SolverUMFPACK::solve: "
              << "SG++ was compiled without UMFPACK support!\n";
    return false;
#endif
}

bool UMFPACK::solve(
        const std::vector<uint32_t> &Ti, const std::vector<uint32_t> &Tj,
        const std::vector<double> &Tx, const std::vector<double> &b,
        std::vector<double> &x, bool initial_output) const
{
#ifdef USEUMFPACK
    if (initial_output)
    {
        tools::printer.printStatusBegin("Solving sparse linear system (UMFPACK)...");
    }
    
    size_t n = b.size();
    
    size_t nnz = Tx.size();
    int *Ti_array = new int[nnz];
    int *Tj_array = new int[nnz];
    double *Tx_array = new double[nnz];
    
    std::copy(Tx.begin(), Tx.end(), Tx_array);
    
    for (size_t k = 0; k < nnz; k++)
    {
        Ti_array[k] = static_cast<int>(Ti[k]);
        Tj_array[k] = static_cast<int>(Tj[k]);
    }
    
    int *Ap = new int[n+1];
    int *Ai = new int[nnz];
    double *Ax = new double[nnz];
    
    tools::printer.printStatusUpdate("step 1: converting to CCS");
    
    if (umfpack_di_triplet_to_col(static_cast<int>(n), static_cast<int>(n), static_cast<int>(nnz),
                                  Ti_array, Tj_array, Tx_array, Ap, Ai, Ax, NULL) != UMFPACK_OK)
    {
        delete[] Ap;
        delete[] Ai;
        delete[] Ax;
        delete[] Ti_array;
        delete[] Tj_array;
        delete[] Tx_array;
        
        tools::printer.printStatusEnd("error: could not convert to CCS!");
        return false;
    }
    
    delete[] Ti_array;
    delete[] Tj_array;
    delete[] Tx_array;
    
    /*if (status_output && (verbosity_level >= 2))
    {
        Output::printStatusUpdate("solving linear system (step 2, " +
                            std::to_string(m) + "x" + std::to_string(n) +
                            ", nnz = " + std::to_string(nnz) + ")");
    }*/
    
    void *symbolic, *numeric;
    double *x_array = new double[n];
    double *b_array = new double[n];
    
    std::copy(b.begin(), b.end(), b_array);
    
    tools::printer.printStatusNewLine();
    tools::printer.printStatusUpdate("step 2: umfpack_di_symbolic");
    
    if (umfpack_di_symbolic(static_cast<int>(n), static_cast<int>(n), Ap, Ai, Ax,
                            &symbolic, NULL, NULL) != UMFPACK_OK)
    {
        delete[] x_array;
        delete[] b_array;
        delete[] Ap;
        delete[] Ai;
        delete[] Ax;
        
        tools::printer.printStatusEnd("error: could solve via umfpack_di_symbolic!");
        return false;
    }
    
    tools::printer.printStatusNewLine();
    tools::printer.printStatusUpdate("step 3: umfpack_di_numeric");
    
    if (umfpack_di_numeric(Ap, Ai, Ax, symbolic, &numeric, NULL, NULL) != UMFPACK_OK)
    {
        umfpack_di_free_symbolic(&symbolic);
        delete[] x_array;
        delete[] b_array;
        delete[] Ap;
        delete[] Ai;
        delete[] Ax;
        
        tools::printer.printStatusEnd("error: could solve via umfpack_di_numeric!");
        return false;
    }
    
    umfpack_di_free_symbolic(&symbolic);
    
    tools::printer.printStatusNewLine();
    tools::printer.printStatusUpdate("step 4: umfpack_di_solve");
    
    if (umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x_array, b_array, numeric, NULL, NULL) !=
        UMFPACK_OK)
    {
        umfpack_di_free_numeric(&numeric);
        delete[] x_array;
        delete[] b_array;
        delete[] Ap;
        delete[] Ai;
        delete[] Ax;
        
        tools::printer.printStatusEnd("error: could solve via umfpack_di_solve!");
        return false;
    }
    
    umfpack_di_free_numeric(&numeric);
    
    x.assign(x_array, x_array + n);
    
    delete[] x_array;
    delete[] b_array;
    delete[] Ap;
    delete[] Ai;
    delete[] Ax;
    
    tools::printer.printStatusEnd();
    
    return true;
#else
    std::cerr << "Error in sg::opt::sle::SolverUMFPACK::solve: "
              << "SG++ was compiled without UMFPACK support!\n";
    return false;
#endif
}

}
}
}
}
