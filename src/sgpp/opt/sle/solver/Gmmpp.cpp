#include "opt/sle/solver/Gmmpp.hpp"
#include "opt/sle/system/Cloneable.hpp"
#include "opt/tools/Printer.hpp"

#ifdef USEGMMPP
#include <gmm/gmm.h>
#endif

#include <cstddef>
#include <iostream>

namespace sg
{
namespace opt
{
namespace sle
{
namespace solver
{

#ifdef USEGMMPP
void callback(const gmm::iteration &iter)
{
    tools::printer.printStatusUpdate("solving with Gmm++ "
            "(k = " + std::to_string(iter.get_iteration()) +
            ", residual norm = " + std::to_string(iter.get_res()) + ")");
}

bool solveInternal(const gmm::csr_matrix<double> &A,
                   const std::vector<double> &b, std::vector<double> &x)
{
    x.assign(b.size(), 0.0);
    gmm::warning_level::level(1);
    
    tools::printer.printStatusUpdate("constructing preconditioner");
    //gmm::identity_matrix P;
    //gmm::diagonal_precond<gmm::csr_matrix<double> > P(A);
    gmm::ilu_precond<gmm::csr_matrix<double> > P(A);
    //gmm::ilut_precond<gmm::csr_matrix<double> > P(A, 10, 1e-4);
    //gmm::ilutp_precond<gmm::csr_matrix<double> > P(A, 10, 1e-4);
    
    tools::printer.printStatusNewLine();
    tools::printer.printStatusUpdate("solving with Gmm++");
    
    gmm::iteration iter(1e-6, 0, 100000);
    iter.set_callback(&callback);
    
    try
    {
        //gmm::bicgstab(A, x, b, P, iter);
        gmm::gmres(A, x, b, P, 50, iter);
        //gmm::qmr(A, x, b, P, iter);
        
        double res = iter.get_res();
        
        if (iter.converged() && (res < 1e3))
        {
            tools::printer.printStatusUpdate("solving with Gmm++ "
                    "(k = " + std::to_string(iter.get_iteration()) +
                    ", residual norm = " + std::to_string(res) + ")");
            tools::printer.printStatusEnd();
            return true;
        } else
        {
            gmm::identity_matrix P;
            
            tools::printer.printStatusNewLine();
            tools::printer.printStatusUpdate(
                    "solving with preconditioner failed, trying again without one");
            tools::printer.printStatusNewLine();
            
            gmm::gmres(A, x, b, P, 50, iter);
            res = iter.get_res();
            
            if (iter.converged() && (res < 1e3))
            {
                tools::printer.printStatusUpdate("solving with Gmm++ "
                        "(k = " + std::to_string(iter.get_iteration()) +
                        ", residual norm = " + std::to_string(res) + ")");
                tools::printer.printStatusEnd();
                return true;
            } else
            {
                tools::printer.printStatusEnd(
                        "error: could not solve linear system, method didn't converge");
                return false;
            }
        }
    } catch (std::exception &e)
    {
        tools::printer.printStatusEnd(
                "error: could not solve linear system, what(): " + std::string(e.what()));
        return false;
    }
}
#endif

bool Gmmpp::solve(system::System &system, std::vector<double> &x) const
{
#ifdef USEGMMPP
    tools::printer.printStatusBegin("Solving linear system (Gmm++)...");
    
    size_t n = system.getDimension();
    const std::vector<double> &b = system.getRHS();
    size_t nnz = 0;
    gmm::csr_matrix<double> A2;
    
    {
        gmm::row_matrix<gmm::rsvector<double> > A(n, n);
        
        std::vector<size_t> Ti;
        std::vector<size_t> Tj;
        std::vector<double> Tx;
        
        #pragma omp parallel if (system.isCloneable()) \
                shared(system, n, Ti, Tj, Tx, A, nnz, tools::printer) default(none)
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
            for (size_t i = 0; i < n; i++)
            {
                for (size_t j = 0; j < n; j++)
                {
                    double entry = system2->getMatrixEntry(i, j);
                    
                    if (entry != 0)
                    {
                        #pragma omp critical
                        {
                            Ti.push_back(i);
                            Tj.push_back(j);
                            Tx.push_back(entry);
                            A(i, j) = entry;
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
        
        gmm::clean(A, 1e-12);
        gmm::copy(A, A2);
        
        // TODO: remove this if not needed anymore
        tools::printer.printVectorToFile("data/Ti.dat", Ti);
        tools::printer.printVectorToFile("data/Tj.dat", Tj);
        tools::printer.printVectorToFile("data/Tx.dat", Tx);
        tools::printer.printVectorToFile("data/b.dat", b);
    }
    
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
    
    bool result = solveInternal(A2, b, x);
    return result;
#else
    std::cerr << "Error in sg::opt::sle::solver::Gmmpp::solve: "
              << "SG++ was compiled without Gmm++ support!\n";
    return false;
#endif
}

}
}
}
}
