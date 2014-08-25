/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_TEST_GRIEWANK_HPP
#define SGPP_OPT_FUNCTION_TEST_GRIEWANK_HPP

#include "opt/function/test/Test.hpp"

#include <cmath>

namespace sg
{
namespace opt
{
namespace function
{
namespace test
{

/**
 * Griewank test function.
 * 
 * Definition:
 * \f$f(\vec{x}) := 1 + \frac{\lVert \vec{x} \rVert_2^2}{4000}
 *                    - \prod_{t=1}^d \cos\!\left(\frac{x_t}{\sqrt{t}}\right)\f$,
 * \f$\vec{x} \in [-600, 600]^d\f$,
 * \f$\vec{x}_{\text{opt}} = \vec{0}\f$,
 * \f$f_{\text{opt}} = 0\f$
 * (domain scaled to \f$[0, 1]^d\f$)
 */
class Griewank : public Test
{
public:
    /**
     * Constructor.
     * 
     * @param d     dimension of the domain
     */
    Griewank(size_t d) : Test(d)
    {
    }
    
    /**
     * Evaluates the test function.
     * 
     * @param x     point \f$\vec{x} \in [0, 1]^d\f$
     * @return      \f$f(\vec{x})\f$
     */
    double evalUndisplaced(const std::vector<double> &x)
    {
        double result = 1.0;
        double tmp = 1.0;
        
        for (size_t t = 0; t < d; t++)
        {
            const double xt = 1200.0 * x[t] - 600.0;
            result += xt*xt / 4000.0;
            tmp *= cos(xt / sqrt(static_cast<double>(t+1)));
        }
        
        result -= tmp;
        return result;
    }
    
    /**
     * Returns minimal point and function value of the test function.
     * 
     * @param[out] x    minimal point \f$\vec{x}_{\text{opt}} \in [0, 1]^d\f$
     * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
     */
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = std::vector<double>(d, 0.5);
        return 0.0;
    }
    
    /**
     * @return clone of the object
     */
    virtual tools::SmartPointer<Objective> clone()
    {
        return tools::SmartPointer<Objective>(new Griewank(*this));
    }
};

}
}
}
}

#endif
