/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/optimizer/NLCG.hpp"
#include "opt/optimizer/LineSearchArmijo.hpp"
#include "opt/tools/Printer.hpp"

#include <numeric>

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
    
    // negated gradient as starting search direction
    for (size_t t = 0; t < d; t++)
    {
        s[t] = -grad_fx[t];
    }
    
    for (k = 0; k < N; k++)
    {
        // exit if norm small enough
        if (grad_fx_norm < tol)
        {
            break;
        }
        
        // normalize search direction
        double s_norm = std::sqrt(std::inner_product(s.begin(), s.end(), s.begin(), 0.0));
        
        for (size_t t = 0; t < d; t++)
        {
            s_normalized[t] = s[t] / s_norm;
        }
        
        // line search
        if (!lineSearchArmijo(*f, beta, gamma, tol, eps, x, fx, grad_fx, s_normalized, y))
        {
            // line search failed ==> exit
            // (either a "real" error occured or the improvement achieved is too small)
            break;
        }
        
        // calculate gradient and norm
        fy = f_gradient->evalGradient(y, grad_fy);
        grad_fy_norm = grad_fy.l2Norm();
        
        double beta = 0.0;
        
        // the method is restarted (beta = 0), if the following criterion is *not* met
        if (std::abs(grad_fy.dotProduct(grad_fx)) / (grad_fy_norm*grad_fy_norm) < alpha)
        {
            // Polak-Ribiere coefficient
            for (size_t t = 0; t < d; t++)
            {
                beta += grad_fy[t] * (grad_fy[t] - grad_fx[t]);
            }
            
            beta /= grad_fx_norm*grad_fx_norm;
        }
        
        // new search direction
        for (size_t t = 0; t < d; t++)
        {
            s[t] = beta * s[t] - grad_fy[t];
        }
        
        // status printing
        tools::printer.printStatusUpdate(toString(k) + " steps, f(x) = " + toString(fx));
        
        x = y;
        fx = fy;
        grad_fx = grad_fy;
        grad_fx_norm = grad_fy_norm;
    }
    
    xopt = x;
    
    tools::printer.printStatusUpdate(toString(k) + " steps, f(x) = " + toString(fx));
    tools::printer.printStatusEnd();
    
    return fx;
}

tools::SmartPointer<Optimizer> NLCG::clone()
{
    tools::SmartPointer<Optimizer> result(
            new NLCG(*f, *f_gradient, N, beta, gamma, tol, eps, alpha));
    result->setStartingPoint(x0);
    return result;
}

const tools::SmartPointer<function::ObjectiveGradient> &NLCG::getObjectiveGradient() const
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
