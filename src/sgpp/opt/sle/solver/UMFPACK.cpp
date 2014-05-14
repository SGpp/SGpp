#include "opt/sle/solver/UMFPACK.hpp"
#include "opt/sle/system/Cloneable.hpp"
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
    
    #pragma omp parallel if (system.isCloneable()) \
            shared(system, n, Ti, Tj, Tx, nnz, tools::printer) default(none)
    {
        system::System *system2;
        std::unique_ptr<system::System> cloned_system;
        
        if (system.isCloneable())
        {
            cloned_system = dynamic_cast<system::Cloneable &>(system).clone();
            system2 = cloned_system.get();
        } else
        {
            system2 = &system;
        }
        
        #pragma omp for ordered schedule(dynamic)
        for (uint32_t i = 0; i < n; i++)
        {
            for (uint32_t j = 0; j < n; j++)
            {
                double entry = system2->getMatrixEntry(i, j);
                
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
    
    typedef UF_long sslong;
    
    size_t n = b.size();
    size_t nnz = Tx.size();
    
    x = std::vector<double>(n, 0.0);
    
    std::unique_ptr<sslong[]> Ap(new sslong[n+1]);
    std::unique_ptr<sslong[]> Ai(new sslong[nnz]);
    std::unique_ptr<double[]> Ax(new double[nnz]);
    
    sslong result;
    
    {
        std::unique_ptr<sslong[]> Ti_array(new sslong[nnz]);
        std::unique_ptr<sslong[]> Tj_array(new sslong[nnz]);
        
        for (size_t k = 0; k < nnz; k++)
        {
            Ti_array[k] = static_cast<sslong>(Ti[k]);
            Tj_array[k] = static_cast<sslong>(Tj[k]);
        }
        
        tools::printer.printStatusUpdate("step 1: umfpack_dl_triplet_to_col");
        
        result = umfpack_dl_triplet_to_col(
                static_cast<sslong>(n), static_cast<sslong>(n), static_cast<sslong>(nnz),
                Ti_array.get(), Tj_array.get(), &Tx[0],
                Ap.get(), Ai.get(), Ax.get(), NULL);
        
        if (result != UMFPACK_OK)
        {
            std::stringstream msg;
            msg << "error: could not convert to CCS via umfpack_dl_triplet_to_col, error code "
                    << result;
            tools::printer.printStatusEnd(msg.str());
            return false;
        }
    }
    
    /*if (status_output && (verbosity_level >= 2))
    {
        Output::printStatusUpdate("solving linear system (step 2, " +
                            std::to_string(m) + "x" + std::to_string(n) +
                            ", nnz = " + std::to_string(nnz) + ")");
    }*/
    
    void *symbolic, *numeric;
    
    tools::printer.printStatusNewLine();
    tools::printer.printStatusUpdate("step 2: umfpack_dl_symbolic");
    
    result = umfpack_dl_symbolic(static_cast<sslong>(n), static_cast<sslong>(n),
                                 Ap.get(), Ai.get(), Ax.get(), &symbolic, NULL, NULL);
    
    if (result != UMFPACK_OK)
    {
        std::stringstream msg;
        msg << "error: could solve via umfpack_dl_symbolic, error code " << result;
        tools::printer.printStatusEnd(msg.str());
        return false;
    }
    
    tools::printer.printStatusNewLine();
    tools::printer.printStatusUpdate("step 3: umfpack_dl_numeric");
    
    result = umfpack_dl_numeric(Ap.get(), Ai.get(), Ax.get(), symbolic,
                                &numeric, NULL, NULL);
    
    if (result != UMFPACK_OK)
    {
        std::stringstream msg;
        msg << "error: could solve via umfpack_dl_numeric, error code " << result;
        tools::printer.printStatusEnd(msg.str());
        umfpack_dl_free_symbolic(&symbolic);
        return false;
    }
    
    umfpack_dl_free_symbolic(&symbolic);
    
    tools::printer.printStatusNewLine();
    tools::printer.printStatusUpdate("step 4: umfpack_dl_solve");
    
    result = umfpack_dl_solve(UMFPACK_A, Ap.get(), Ai.get(), Ax.get(), &x[0], &b[0],
                              numeric, NULL, NULL);
    
    if (result != UMFPACK_OK)
    {
        std::stringstream msg;
        msg << "error: could solve via umfpack_dl_solve, error code " << result;
        tools::printer.printStatusEnd(msg.str());
        umfpack_dl_free_numeric(&numeric);
        return false;
    }
    
    umfpack_dl_free_numeric(&numeric);
    tools::printer.printStatusEnd();
    
    return true;
#else
    std::cerr << "Error in sg::opt::sle::solver::UMFPACK::solve: "
              << "SG++ was compiled without UMFPACK support!\n";
    return false;
#endif
}

}
}
}
}
