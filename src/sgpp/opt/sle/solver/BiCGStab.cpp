#include "opt/sle/solver/BiCGStab.hpp"
#include "opt/tools/Printer.hpp"

#include <cmath>
#include <numeric>

namespace sg
{
namespace opt
{
namespace sle
{
namespace solver
{

const double BiCGStab::DEFAULT_TOLERANCE = 1e-10;

bool solveInternal(system::System &system, size_t N, double tol,
                   const std::vector<double> &x0, std::vector<double> &x)
{
    size_t n = system.getDimension();
    const std::vector<double> &b = system.getRHS();
    std::vector<double> r(n, 0.0);
    //double b_norm_squared = std::inner_product(b.begin(), b.end(), b.begin(), 0.0);
    
    std::vector<double> my_x0 = x0;
    
    if (my_x0.empty())
    {
        my_x0 = std::vector<double>(n, 0.0);
    }
    
    x = std::vector<double>(n, 0.0);
    system.matrixVectorMultiplication(my_x0, r);
    
    for (size_t i = 0; i < n; i++)
    {
        r[i] = b[i] - r[i];
    }
    
    std::vector<double> r0hat(r);
    double rho = 1.0;
    double alpha = 1.0;
    double omega = 1.0;
    std::vector<double> v(n, 0.0);
    std::vector<double> p(n, 0.0);
    std::vector<double> s(n, 0.0);
    std::vector<double> t(n, 0.0);
    double r_norm_squared = 0.0;
    
    for (size_t k = 0; k < N; k++)
    {
        double last_rho = rho;
        rho = std::inner_product(r0hat.begin(), r0hat.end(), r.begin(), 0.0);
        double beta = (rho / last_rho) * (alpha / omega);
        
        for (size_t i = 0; i < n; i++)
        {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        
        system.matrixVectorMultiplication(p, v);
        alpha = rho / std::inner_product(r0hat.begin(), r0hat.end(), v.begin(), 0.0);
        
        for (size_t i = 0; i < n; i++)
        {
            s[i] = r[i] - alpha * v[i];
        }
        
        system.matrixVectorMultiplication(s, t);
        omega = std::inner_product(t.begin(), t.end(), s.begin(), 0.0) /
                std::inner_product(t.begin(), t.end(), t.begin(), 0.0);
        //if (k <= 10) std::cout << "\n" << k << " " << std::inner_product(t.begin(), t.end(), t.begin(), 0.0) << "\n";
        
        if (std::isnan(omega))
        {
            return false;
        }
        
        for (size_t i = 0; i < n; i++)
        {
            x[i] = x[i] + alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * t[i];
        }
        
        r_norm_squared = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
        
        /*if (status_output && (verbosity_level >= 2) && (k % 10 == 0))
        {
            std::stringstream msg;
            msg << "residual norm = " << sqrt(r_norm_squared);
            Output::printStatusUpdate(msg.str());
        }*/
        tools::printer.printStatusUpdate("k = " + std::to_string(k) + ", residual norm = " +
                                         std::to_string(sqrt(r_norm_squared)));
        
        //if (r_norm_squared / b_norm_squared < tol*tol)
        if (r_norm_squared < tol*tol)
        {
            break;
        }
    }
    
    tools::printer.printStatusUpdate("residual norm = " + std::to_string(sqrt(r_norm_squared)));
    
    return true;
}

BiCGStab::BiCGStab() : BiCGStab(DEFAULT_MAX_IT_COUNT, DEFAULT_TOLERANCE, std::vector<double>())
{
}

BiCGStab::BiCGStab(size_t max_it_count, double tolerance, const std::vector<double> &x0) :
    Solver(),
    N(max_it_count),
    tol(tolerance),
    x0(x0)
{
}

bool BiCGStab::solve(system::System &system, std::vector<double> &x) const
{
    tools::printer.printStatusBegin("Solving linear system (BiCGStab)...");
    
    bool result = solveInternal(system, N, tol, x0, x);
    
    if (result)
    {
        tools::printer.printStatusEnd();
        return true;
    } else
    {
        tools::printer.printStatusEnd("error: could not solve linear system!");
        return false;
    }
}

/*bool SolverBiCGStab::solve(const std::vector<uint32_t> &Ti, const std::vector<uint32_t> &Tj,
                           const std::vector<double> &Tx, const std::vector<double> &b,
                           std::vector<double> &x) const
{
    tools::printer.printStatusBegin("Solving sparse linear system (BiCGStab)...");
    
    size_t nnz = Tx.size();
    size_t n = b.size();
    base::DataMatrix A(n, n);
    
    A.setAll(0.0);
    
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
        
        A.set(Ti[k], Tj[k], Tx[k]);
    }
    
    tools::printer.printStatusUpdate("constructing full matrix (100.0%)");
    tools::printer.printStatusNewLine();
    
    std::stringstream msg;
    double nnz_ratio = static_cast<double>(nnz) /
                       (static_cast<double>(n) * static_cast<double>(n));
    msg << "nnz ratio: " << static_cast<double>(static_cast<int>(nnz_ratio * 1000.0)) / 10.0
        << "%";
    tools::printer.printStatusUpdate(msg.str());
    tools::printer.printStatusNewLine();
    
    FullSystem system(n, A, b);
    
    bool result = solveInternal(system, N, tol, x0, x);
    
    if (result)
    {
        tools::printer.printStatusEnd();
        return true;
    } else
    {
        tools::printer.printStatusEnd("error: could not solve linear system!");
        return false;
    }
}*/

size_t BiCGStab::getMaxItCount() const
{
    return N;
}

void BiCGStab::setMaxItCount(size_t max_it_count)
{
    N = max_it_count;
}

double BiCGStab::getTolerance() const
{
    return tol;
}

void BiCGStab::setTolerance(double tolerance)
{
    tol = tolerance;
}

const std::vector<double> &BiCGStab::getX0() const
{
    return x0;
}

void BiCGStab::setX0(const std::vector<double> &x0)
{
    this->x0 = x0;
}

}
}
}
}
