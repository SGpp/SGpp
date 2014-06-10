#include "opt/optimizer/GradientMethod.hpp"
//#include "opt/optimizer/LineSearchGolden.hpp"
#include "opt/optimizer/LineSearchArmijo.hpp"
//#include "opt/optimizer/LineSearchPowellWolfe.hpp"
//#include "opt/optimizer/LineSearchNewton.hpp"
#include "opt/tools/Printer.hpp"

namespace sg
{
namespace opt
{
namespace optimizer
{

const double GradientMethod::DEFAULT_BETA = 0.5;
const double GradientMethod::DEFAULT_GAMMA = 1e-2;
const double GradientMethod::DEFAULT_TOLERANCE = 1e-8;
const double GradientMethod::DEFAULT_EPSILON = 1e-18;

GradientMethod::GradientMethod(
        function::Objective &f,
        function::ObjectiveGradient &f_gradient,
        size_t N, double beta, double gamma, double tolerance, double epsilon) :
    Optimizer(f, N),
    f_gradient(f_gradient.clone()),
    beta(beta),
    gamma(gamma),
    tol(tolerance),
    eps(epsilon)
{
}

double GradientMethod::optimize(std::vector<double> &xopt)
{
    tools::printer.printStatusBegin("Optimizing (gradient method)...");
    
    size_t d = f->getDimension();
    std::vector<double> x(x0);
    double fx = 0.0;
    
    base::DataVector grad_fx(d);
    std::vector<double> s(d, 0.0);
    std::vector<double> y(d, 0.0);
    size_t k;
    
    for (k = 0; k < N; k++)
    {
        fx = f_gradient->evalGradient(x, grad_fx);
        
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
        {
            std::stringstream msg;
            msg << k << " steps, f(x) = " << fx;
            //msg << k << " steps, f(x) = " << fx << ", |grad f(x)| = " << grad_fx_norm;
            tools::printer.printStatusUpdate(msg.str());
        }
        
        //if (!lineSearchGolden(*f, tol, x, fx, s, y))
        if (!lineSearchArmijo(*f, beta, gamma, tol, eps, x, fx, grad_fx, s, y))
        //if (!lineSearchPowellWolfe(*f, *f_gradient, 1e-4, 0.1, x, fx, grad_fx, s, y))
        //if (!lineSearchNewton(*f_hessian, tol, x, s, y))
        {
            break;
        }
        
        x = y;
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

std::unique_ptr<Optimizer> GradientMethod::clone()
{
    std::unique_ptr<Optimizer> result(
            new GradientMethod(*f, *f_gradient, N, beta, gamma, tol, eps));
    result->setStartingPoint(x0);
    return result;
}

const std::unique_ptr<function::ObjectiveGradient> &GradientMethod::getObjectiveGradient() const
{
    return f_gradient;
}

double GradientMethod::getBeta() const
{
    return beta;
}

void GradientMethod::setBeta(double beta)
{
    this->beta = beta;
}

double GradientMethod::getGamma() const
{
    return gamma;
}

void GradientMethod::setGamma(double gamma)
{
    this->gamma = gamma;
}

double GradientMethod::getTolerance() const
{
    return tol;
}

void GradientMethod::setTolerance(double tolerance)
{
    tol = tolerance;
}

double GradientMethod::getEpsilon() const
{
    return eps;
}

void GradientMethod::setEpsilon(double epsilon)
{
    eps = epsilon;
}

}
}
}
