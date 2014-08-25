/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_OPTIMIZER_LINESEARCHARMIJO_HPP
#define SGPP_OPT_OPTIMIZER_LINESEARCHARMIJO_HPP

#include "base/datatypes/DataVector.hpp"
#include "opt/function/Objective.hpp"

#include <cstddef>
#include <cmath>

namespace sg
{
namespace opt
{
namespace optimizer
{

/**
 * Line search (1D optimization on a line) with Armijo's rule used in gradient-based optimization.
 * 
 * Armijo's rule calculates \f$\sigma = \beta^k\f$ for \f$k = 0, 1, \dotsc\f$
 * for a fixed \f$\beta \in (0, 1)\f$ and checks if \f$\vec{y} = \vec{x} + \sigma\vec{s}\f$
 * lies in \f$[0, 1]^d\f$ and whether
 * the objective function value improvement meets the condition
 * \f$f(\vec{x}) - f(\vec{y}) \ge \gamma\sigma (-\nabla f(\vec{x}) \cdot \vec{s})\f$
 * for \f$\gamma \in (0, 1)\f$ fixed.
 * 
 * The return value states whether the relative improvement (depending on two tolerances)
 * is big enough to continue the optimization algorithm.
 * 
 * @param       f       objective function
 * @param       beta    \f$\beta \in (0, 1)\f$
 * @param       gamma   \f$\gamma \in (0, 1)\f$
 * @param       tol     tolerance 1 (positive)
 * @param       eps     tolerance 2 (positive)
 * @param       x       point to start the line search in
 * @param       fx      objective function value in x
 * @param       grad_fx objective function gradient in x
 * @param       s       search direction (should be normalized)
 * @param[out]  y       new point, must have the same size as x before calling this function
 * @return              whether the new point reaches an acceptable improvement
 */
inline bool lineSearchArmijo(function::Objective &f, double beta, double gamma,
                             double tol, double eps,
                             const std::vector<double> &x, double fx, base::DataVector &grad_fx,
                             const std::vector<double> &s, std::vector<double> &y)
{
    const size_t d = x.size();
    double sigma = 1.0;
    double ip = 0.0;
    double fy = fx;
    
    // inner product between grad_fx and s
    for (size_t t = 0; t < d; t++)
    {
        ip += grad_fx[t] * s[t];
    }
    
    for (size_t k = 0; k < 100; k++)
    {
        bool y_in_domain = true;
        
        // calculate new point
        for (size_t t = 0; t < d; t++)
        {
            y[t] = x[t] + sigma * s[t];
            
            if ((y[t] < 0.0) || (y[t] > 1.0))
            {
                y_in_domain = false;
                break;
            }
        }
        
        // check if y lies in [0, 1]^d
        if (y_in_domain)
        {
            fy = f.eval(y);
            
            double improvement = fx - fy;
            double rhs = gamma * sigma * (-ip);
            
            // check if the absolute improvement is big enough
            if (improvement >= rhs)
            {
                return (std::abs(fx - fy) >= tol * (std::abs(fx) + std::abs(fy) + eps));
            }
        }
        
        // next sigma
        sigma *= beta;
    }
    
    return false;
}

}
}
}

#endif
