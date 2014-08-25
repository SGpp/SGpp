/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

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
/**
 * Gmm++ status printing callback.
 * 
 * @param iter  iteration information
 */
void callback(const gmm::iteration &iter)
{
    tools::printer.printStatusUpdate("solving with Gmm++ "
            "(k = " + toString(iter.get_iteration()) +
            ", residual norm = " + toString(iter.get_res()) + ")");
}

/**
 * @param       A   coefficient matrix
 * @param       b   right-hand side
 * @param[out]  x   solution of the linear system
 * @return          whether all went well (false if errors occurred)
 */
bool solveInternal(const gmm::csr_matrix<double> &A,
                   const std::vector<double> &b, std::vector<double> &x)
{
    // allow warnings
    gmm::warning_level::level(1);
    x.assign(b.size(), 0.0);
    
    // ILU preconditioning
    tools::printer.printStatusUpdate("constructing preconditioner");
    gmm::ilu_precond<gmm::csr_matrix<double> > P(A);
    
    tools::printer.printStatusNewLine();
    tools::printer.printStatusUpdate("solving with Gmm++");
    
    gmm::iteration iter(1e-6, 0, 100000);
    iter.set_callback(&callback);
    
    try
    {
        // call GMRES
        gmm::gmres(A, x, b, P, 50, iter);
        
        double res = iter.get_res();
        
        if (iter.converged() && (res < 1e3))
        {
            // GMRES converged
            tools::printer.printStatusUpdate("solving with Gmm++ "
                    "(k = " + toString(iter.get_iteration()) +
                    ", residual norm = " + toString(res) + ")");
            tools::printer.printStatusEnd();
            return true;
        } else
        {
            // GMRES didn't converge ==> try again without preconditioner
            gmm::identity_matrix P;
            
            tools::printer.printStatusNewLine();
            tools::printer.printStatusUpdate(
                    "solving with preconditioner failed, trying again without one");
            tools::printer.printStatusNewLine();
            
            // call GMRES again
            gmm::gmres(A, x, b, P, 50, iter);
            res = iter.get_res();
            
            if (iter.converged() && (res < 1e3))
            {
                tools::printer.printStatusUpdate("solving with Gmm++ "
                        "(k = " + toString(iter.get_iteration()) +
                        ", residual norm = " + toString(res) + ")");
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

bool Gmmpp::solve(system::System &system, const std::vector<double> &b,
                  std::vector<double> &x) const
{
#ifdef USEGMMPP
    tools::printer.printStatusBegin("Solving linear system (Gmm++)...");
    
    const size_t n = system.getDimension();
    size_t nnz = 0;
    gmm::csr_matrix<double> A2;
    
    {
        gmm::row_matrix<gmm::rsvector<double> > A(n, n);
        
        // parallelize only if the system is cloneable
        #pragma omp parallel if (system.isCloneable()) \
                shared(system, A, nnz, tools::printer) default(none)
        {
            system::System *system2;
            tools::SmartPointer<system::System> cloned_system;
            
            if (system.isCloneable())
            {
                cloned_system = dynamic_cast<system::Cloneable &>(system).clone();
                system2 = cloned_system.get();
            } else
            {
                system2 = &system;
            }
            
            // copy system matrix to Gmm++ matrix object
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
                            A(i, j) = entry;
                            nnz++;
                        }
                    }
                }
                
                // status message
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
        
        // the Gmm++ manual said to do so... no idea if that's necessary
        // (but I'm a big fan of RTFM ;) )
        gmm::clean(A, 1e-12);
        gmm::copy(A, A2);
    }
    
    tools::printer.printStatusUpdate("constructing sparse matrix (100.0%)");
    tools::printer.printStatusNewLine();
    
    // print ratio of nonzero entries
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
