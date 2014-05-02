#include "opt/optimization/OptimizerGradientMethod.hpp"
#include "opt/optimization/ArmijoRule.hpp"
#include "opt/tools/Printer.hpp"

namespace sg
{
namespace opt
{
namespace optimization
{

const double OptimizerGradientMethod::DEFAULT_BETA = 0.5;
const double OptimizerGradientMethod::DEFAULT_GAMMA = 1e-2;
const double OptimizerGradientMethod::DEFAULT_TOLERANCE = 1e-10;

OptimizerGradientMethod::OptimizerGradientMethod(
        function::ObjectiveFunction &f,
        function::ObjectiveFunctionGradient &f_gradient) :
    OptimizerGradientMethod(f, f_gradient, DEFAULT_MAX_IT_COUNT,
                            DEFAULT_BETA, DEFAULT_GAMMA, DEFAULT_TOLERANCE)
{
}

OptimizerGradientMethod::OptimizerGradientMethod(
        function::ObjectiveFunction &f,
        function::ObjectiveFunctionGradient &f_gradient,
        size_t N, double beta, double gamma, double tolerance) :
    Optimizer(f, N),
    f_gradient(f_gradient),
    beta(beta),
    gamma(gamma),
    tol(tolerance)
{
}

void OptimizerGradientMethod::optimize(std::vector<double> &xopt)
{
    tools::printer.printStatusBegin("Optimizing (gradient method)...");
    
    size_t d = f.getDimension();
    std::vector<double> x(x0);
    double fx;
    
    base::DataVector grad_fx(d);
    std::vector<double> s(d, 0.0);
    std::vector<double> y(d, 0.0);
    size_t k;
    
    for (k = 0; k < N; k++)
    {
        fx = f_gradient.evalGradient(x, grad_fx);
        
        double grad_fx_norm = grad_fx.l2Norm();
        
        if (grad_fx_norm < tol)
        {
            break;
        }
        
        for (size_t t = 0; t < d; t++)
        {
            s[t] = -grad_fx[t] / grad_fx_norm;
        }
        
        /*std::cout << "\nk: " << k << "\n";
        std::cout << "x: [" << x[0] << ", " << x[1] << "]\n";
        std::cout << "fx: " << fx << "\n";
        std::cout << "grad_fx: [" << grad_fx[0] << ", " << grad_fx[1] << "]\n";
        std::cout << "s: [" << s[0] << ", " << s[1] << "]\n";*/
        if (k % 10 == 0)
        {
            std::stringstream msg;
            msg << k << " steps, f(x) = " << fx;
            //msg << k << " steps, f(x) = " << fx << ", |grad f(x)| = " << grad_fx_norm;
            tools::printer.printStatusUpdate(msg.str());
        }
        
        if (!armijoRule(f, beta, gamma, x, fx, grad_fx, s, y))
        {
            break;
        }
        
        x = y;
    }
    
    xopt = x;
    
    std::stringstream msg;
    msg << k << " steps, f(x) = " << fx;
    tools::printer.printStatusUpdate(msg.str());
    tools::printer.printStatusEnd();
}

function::ObjectiveFunctionGradient
        &OptimizerGradientMethod::getObjectiveFunctionGradient() const
{
    return f_gradient;
}

double OptimizerGradientMethod::getBeta() const
{
    return beta;
}

void OptimizerGradientMethod::setBeta(double beta)
{
    this->beta = beta;
}

double OptimizerGradientMethod::getGamma() const
{
    return gamma;
}

void OptimizerGradientMethod::setGamma(double gamma)
{
    this->gamma = gamma;
}

double OptimizerGradientMethod::getTolerance() const
{
    return tol;
}

void OptimizerGradientMethod::setTolerance(double tolerance)
{
    tol = tolerance;
}

}
}
}
