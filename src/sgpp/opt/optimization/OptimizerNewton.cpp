#include "opt/optimization/OptimizerNewton.hpp"
#include "opt/optimization/ArmijoRule.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "opt/sle/system/Full.hpp"
#include "opt/tools/Printer.hpp"

#include <algorithm>
#include <numeric>

namespace sg
{
namespace opt
{
namespace optimization
{

const double OptimizerNewton::DEFAULT_ALPHA1 = 1e-6;
const double OptimizerNewton::DEFAULT_ALPHA2 = 1e-6;
const double OptimizerNewton::DEFAULT_BETA = 0.5;
const double OptimizerNewton::DEFAULT_GAMMA = 1e-2;
const double OptimizerNewton::DEFAULT_P = 0.1;
const double OptimizerNewton::DEFAULT_TOLERANCE = 1e-10;

OptimizerNewton::OptimizerNewton(
        function::ObjectiveFunction &f,
        function::ObjectiveFunctionHessian &f_hessian) :
    OptimizerNewton(f, f_hessian, DEFAULT_MAX_IT_COUNT,
                    DEFAULT_ALPHA1, DEFAULT_ALPHA2, DEFAULT_BETA, DEFAULT_GAMMA,
                    DEFAULT_P, DEFAULT_TOLERANCE, default_sle_solver)
{
}

OptimizerNewton::OptimizerNewton(
        function::ObjectiveFunction &f,
        function::ObjectiveFunctionHessian &f_hessian,
        size_t max_it_count, double alpha1, double alpha2, double beta, double gamma,
        double p, double tolerance, const sle::solver::Solver &sle_solver) :
    Optimizer(f, max_it_count),
    f_hessian(f_hessian),
    alpha1(alpha1),
    alpha2(alpha2),
    beta(beta),
    gamma(gamma),
    p(p),
    tol(tolerance),
    default_sle_solver(sle::solver::BiCGStab()),
    sle_solver(sle_solver)
{
}

void OptimizerNewton::optimize(std::vector<double> &xopt)
{
    tools::printer.printStatusBegin("Optimizing (Newton)...");
    
    size_t d = f.getDimension();
    std::vector<double> x = x0;
    double fx;
    
    std::vector<double> ls_solver_x0(d, 0.0);
    bool ls_solved;
    std::vector<double> dk(d, 0.0);
    
    base::DataVector grad_fx(d);
    base::DataMatrix hessian_fx(d, d);
    std::vector<double> s(d, 0.0);
    std::vector<double> y(d, 0.0);
    size_t verbosity = tools::printer.getVerbosity();
    
    sle::system::Full system(hessian_fx, s);
    size_t k;
    
    for (k = 0; k < N; k++)
    {
        fx = f_hessian.evalHessian(x, grad_fx, hessian_fx);
        
        double grad_fx_norm = grad_fx.l2Norm();
        
        if (grad_fx_norm < tol)
        {
            break;
        }
        
        for (size_t t = 0; t < d; t++)
        {
            s[t] = -grad_fx[t];
        }
        
        /*ls_solved = Coeffs::BiCGSTAB::solveLinearSystem(hessian_fx, s,
                                                        ls_solver_x0, ls_solver_tol, false, dk);
        double dk_norm = sqrt(std::inner_product(dk.begin(), dk.end(), dk.begin(), 0.0));*/
        
        system.setA(hessian_fx);
        system.setRHS(s);
        tools::printer.setVerbosity(0);
        ls_solved = sle_solver.solve(system, dk);
        tools::printer.setVerbosity(verbosity);
        
        double dk_norm = sqrt(std::inner_product(dk.begin(), dk.end(), dk.begin(), 0.0));
        
        if (ls_solved && (std::inner_product(s.begin(), s.end(), dk.begin(), 0.0) >=
                          std::min(alpha1, alpha2*std::pow(dk_norm, p)) * dk_norm*dk_norm))
        {
            for (size_t t = 0; t < d; t++)
            {
                s[t] = dk[t] / dk_norm;
            }
        } else
        {
            for (size_t t = 0; t < d; t++)
            {
                s[t] = s[t] / grad_fx_norm;
            }
        }
        
        if (k % 10 == 0)
        {
            std::stringstream msg;
            msg << k << " steps, f(x) = " << fx;
            //msg << k << " steps, f(x) = " << fx << ", |grad f(x)| = " << grad_fx_norm;
            tools::printer.printStatusUpdate(msg.str());
            //std::cout << "\n" << x << ", " << fx << "\n";
        }
        
        if (!armijoRule(f, beta, gamma, x, fx, grad_fx, s, y))
        {
            //std::cout << "\nBOOM\n";
            break;
        }
        
        x = y;
    }
    
    xopt = x;
    tools::printer.setVerbosity(verbosity);
    
    std::stringstream msg;
    msg << k << " steps, f(x) = " << fx;
    tools::printer.printStatusUpdate(msg.str());
    tools::printer.printStatusEnd();
}

function::ObjectiveFunctionHessian &OptimizerNewton::getObjectiveFunctionHessian() const
{
    return f_hessian;
}

double OptimizerNewton::getAlpha1() const
{
    return alpha1;
}

void OptimizerNewton::setAlpha1(double alpha1)
{
    this->alpha1 = alpha1;
}

double OptimizerNewton::getAlpha2() const
{
    return alpha2;
}

void OptimizerNewton::setAlpha2(double alpha2)
{
    this->alpha2 = alpha2;
}

double OptimizerNewton::getBeta() const
{
    return beta;
}

void OptimizerNewton::setBeta(double beta)
{
    this->beta = beta;
}

double OptimizerNewton::getGamma() const
{
    return gamma;
}

void OptimizerNewton::setGamma(double gamma)
{
    this->gamma = gamma;
}

double OptimizerNewton::getP() const
{
    return p;
}

void OptimizerNewton::setP(double p)
{
    this->p = p;
}

double OptimizerNewton::getTolerance() const
{
    return tol;
}

void OptimizerNewton::setTolerance(double tolerance)
{
    tol = tolerance;
}

}
}
}
