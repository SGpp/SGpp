#include "opt/sle/solver/Eigen.hpp"
#include "opt/sle/system/Cloneable.hpp"
#include "opt/tools/Printer.hpp"

#ifdef USEEIGEN
#include <eigen3/Eigen/Dense>
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

#ifdef USEEIGEN
bool solveInternal(const ::Eigen::MatrixXd &A, const std::vector<double> &b,
                   std::vector<double> &x)
{
    tools::printer.printStatusUpdate("step 1: Householder QR factorization");
    
    ::Eigen::HouseholderQR< ::Eigen::MatrixXd> A_QR = A.householderQr();
    //::Eigen::FullPivLU< ::Eigen::MatrixXd> A_LU = A.fullPivLu();
    
    tools::printer.printStatusNewLine();
    tools::printer.printStatusUpdate("step 2: solving");
    
    ::Eigen::VectorXd b_eigen = ::Eigen::VectorXd::Map(&b[0], b.size());
    ::Eigen::VectorXd x_eigen = A_QR.solve(b_eigen);
    //Eigen::VectorXd x_eigen = A_LU.solve(b_eigen);
    
    if ((A*x_eigen).isApprox(b_eigen))
    {
        x = std::vector<double>(x_eigen.data(), x_eigen.data() + x_eigen.size());
        tools::printer.printStatusEnd();
        return true;
    } else
    {
        tools::printer.printStatusEnd("error: could not solve linear system!");
        return false;
    }
}
#endif

bool Eigen::solve(system::System &system, std::vector<double> &x) const
{
#ifdef USEEIGEN
    tools::printer.printStatusBegin("Solving linear system (Eigen)...");
    
    size_t n = system.getDimension();
    const std::vector<double> &b = system.getRHS();
    ::Eigen::MatrixXd A = ::Eigen::MatrixXd::Zero(n, n);
    size_t nnz = 0;
    
    #pragma omp parallel if (system.isCloneable()) \
            shared(system, n, A, nnz, tools::printer) default(none)
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
                A(i,j) = system2->getMatrixEntry(i, j);
                
                if (A(i,j) != 0)
                {
                    #pragma omp atomic
                    nnz++;
                }
            }
            
            if (i % 100 == 0)
            {
                #pragma omp ordered
                {
                    char str[10];
                    snprintf(str, 10, "%.1f%%",
                             static_cast<double>(i) / static_cast<double>(n) * 100.0);
                    tools::printer.printStatusUpdate("constructing matrix (" +
                                                     std::string(str) + ")");
                }
            }
        }
    }
    
    tools::printer.printStatusUpdate("constructing matrix (100.0%)");
    tools::printer.printStatusNewLine();
    
    {
        char str[10];
        double nnz_ratio = static_cast<double>(nnz) /
                           (static_cast<double>(n) * static_cast<double>(n));
        snprintf(str, 10, "%.1f%%", nnz_ratio * 100.0);
        tools::printer.printStatusUpdate("nnz ratio: " + std::string(str));
        tools::printer.printStatusNewLine();
    }
    
    bool result = solveInternal(A, b, x);
    return result;
#else
    std::cerr << "Error in sg::opt::sle::solver::Eigen::solve: "
              << "SG++ was compiled without Eigen support!\n";
    return false;
#endif
}

/*bool SolverEigen::solve(const std::vector<uint32_t> &Ti, const std::vector<uint32_t> &Tj,
                        const std::vector<double> &Tx, const std::vector<double> &b,
                        std::vector<double> &x) const
{
#ifdef USEEIGEN
    tools::printer.printStatusBegin("Solving linear system (Eigen)...");
    
    size_t nnz = Tx.size();
    size_t n = b.size();
    ::Eigen::MatrixXd A = ::Eigen::MatrixXd::Zero(n, n);
    
    for (size_t k = 0; k < nnz; k++)
    {
        if (k % 100 == 0)
        {
            char str[10];
            snprintf(str, 10, "%.1f%%",
                     static_cast<double>(k) / static_cast<double>(nnz) * 100.0);
            tools::printer.printStatusUpdate("constructing full matrix (" +
                                             std::string(str) + ")");
        }
        
        A(Ti[k], Tj[k]) = Tx[k];
    }
    
    tools::printer.printStatusUpdate("constructing full matrix (100.0%)");
    tools::printer.printStatusNewLine();
    
    bool result = solveInternal(A, b, x);
    
    if (result)
    {
        tools::printer.printStatusEnd();
        return true;
    } else
    {
        tools::printer.printStatusEnd("error: could not solve linear system!");
        return false;
    }
#else
    std::cerr << "Error in sg::opt::sle::SolverEigen::solve: "
              << "SG++ was compiled without Eigen support!\n";
    return false;
#endif
}*/

}
}
}
}
