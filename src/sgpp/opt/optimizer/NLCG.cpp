#include "opt/optimizer/NLCG.hpp"
//#include "opt/optimizer/LineSearchGolden.hpp"
#include "opt/optimizer/LineSearchArmijo.hpp"
//#include "opt/optimizer/LineSearchNewton.hpp"
//#include "opt/optimizer/LineSearchPowellWolfe.hpp"
#include "opt/tools/Printer.hpp"

namespace sg
{
namespace opt
{
namespace optimizer
{

const double NLCG::DEFAULT_BETA = 0.5;
const double NLCG::DEFAULT_GAMMA = 1e-2;
const double NLCG::DEFAULT_TOLERANCE = 1e-8;
const double NLCG::DEFAULT_EPSILON = 1e-18;
const double NLCG::DEFAULT_RESTART_THRESHOLD = 0.1;

NLCG::NLCG(function::Objective &f, function::ObjectiveGradient &f_gradient,
           size_t max_it_count, double beta, double gamma,
           double tolerance, double epsilon, double restart_threshold) :
    Optimizer(f, max_it_count),
    f_gradient(f_gradient.clone()),
    beta(beta),
    gamma(gamma),
    tol(tolerance),
    eps(epsilon),
    alpha(restart_threshold)
{
}

double NLCG::optimize(std::vector<double> &xopt)
{
    tools::printer.printStatusBegin("Optimizing (NLCG)...");
    
    size_t d = f->getDimension();
    std::vector<double> x(x0);
    double fx;
    double fy;
    
    base::DataVector grad_fx(d);
    base::DataVector grad_fy(d);
    std::vector<double> s(d, 0.0);
    std::vector<double> s_normalized(d, 0.0);
    std::vector<double> y(d, 0.0);
    size_t k;
    
    fx = f_gradient->evalGradient(x0, grad_fx);
    double grad_fx_norm = grad_fx.l2Norm();
    double grad_fy_norm = 0.0;
    
    for (size_t t = 0; t < d; t++)
    {
        s[t] = -grad_fx[t];
    }
    
    for (k = 0; k < N; k++)
    {
        if (grad_fx_norm < tol)
        {
            break;
        }
        
        double s_norm = std::sqrt(std::inner_product(s.begin(), s.end(), s.begin(), 0.0));
        
        for (size_t t = 0; t < d; t++)
        {
            s_normalized[t] = s[t] / s_norm;
        }
        
        //if (!lineSearchGolden(*f, tol, x, fx, s, y))
        if (!lineSearchArmijo(*f, beta, gamma, tol, eps, x, fx, grad_fx, s_normalized, y))
        //if (!lineSearchPowellWolfe(*f, *f_gradient, 1e-4, 0.1, x, fx, grad_fx, s_normalized, y))
        //if (!lineSearchNewton(*f_hessian, tol, x, s_normalized, y))
        {
            break;
        }
        
        fy = f_gradient->evalGradient(y, grad_fy);
        grad_fy_norm = grad_fy.l2Norm();
        
        /*if (std::abs(fx - fy) <= 1e-8 * (std::abs(fx) + std::abs(fy) + 1e-18))
        {
            break;
        }*/
        
        double beta = 0.0;
        
        if (std::abs(grad_fy.dotProduct(grad_fx)) / (grad_fy_norm*grad_fy_norm) < alpha)
        {
            for (size_t t = 0; t < d; t++)
            {
                beta += grad_fy[t] * (grad_fy[t] - grad_fx[t]);
            }
            
            beta /= grad_fx_norm*grad_fx_norm;
        }
        
        for (size_t t = 0; t < d; t++)
        {
            s[t] = beta * s[t] - grad_fy[t];
        }
        
        {
            std::stringstream msg;
            msg << k << " steps, f(x) = " << fx;
            //msg << k << " steps, f(x) = " << fx << ", |grad f(x)| = " << grad_fx_norm;
            tools::printer.printStatusUpdate(msg.str());
        }
        
        x = y;
        fx = fy;
        grad_fx = grad_fy;
        grad_fx_norm = grad_fy_norm;
    }
    
    xopt = x;
    
    {
        std::stringstream msg;
        msg << k << " steps, f(x) = " << fx;
        tools::printer.printStatusUpdate(msg.str());
        tools::printer.printStatusEnd();
    }
    
    return fx;
}

std::unique_ptr<Optimizer> NLCG::clone()
{
    std::unique_ptr<Optimizer> result(new NLCG(*f, *f_gradient, N, tol, eps));
    result->setStartingPoint(x0);
    return result;
}

const std::unique_ptr<function::ObjectiveGradient> &NLCG::getObjectiveGradient() const
{
    return f_gradient;
}

double NLCG::getBeta() const
{
    return beta;
}

void NLCG::setBeta(double beta)
{
    this->beta = beta;
}

double NLCG::getGamma() const
{
    return gamma;
}

void NLCG::setGamma(double gamma)
{
    this->gamma = gamma;
}

double NLCG::getTolerance() const
{
    return tol;
}

void NLCG::setTolerance(double tolerance)
{
    tol = tolerance;
}

double NLCG::getEpsilon() const
{
    return eps;
}

void NLCG::setEpsilon(double epsilon)
{
    eps = epsilon;
}

double NLCG::getRestartThreshold() const
{
    return alpha;
}

void NLCG::setRestartThreshold(double restart_threshold)
{
    alpha = restart_threshold;
}

}
}
}
