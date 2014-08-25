/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_TEST_GOLDSTEINPRICE_HPP
#define SGPP_OPT_FUNCTION_TEST_GOLDSTEINPRICE_HPP

#include "opt/function/test/Test.hpp"

namespace sg
{
namespace opt
{
namespace function
{
namespace test
{

/**
 * Goldstein-Price test function.
 * 
 * Definition:
 * \f$f(\vec{x}) := (1 + (x_1+x_2+1)^2 (19 - 14x_1 + 3x_1^2 - 14x_2 + 6x_1 x_2 + 3x_2^2))
 *            \cdot (30 + (2x_1 - 3x_2)^2 (18 - 32x_1 + 12x_1^2 + 48x_2 - 36x_1 x_2 + 27x_2^2))\f$,
 * \f$\vec{x} \in [-2, 2]^2\f$,
 * \f$\vec{x}_{\text{opt}} = (0, -1)^{\mathrm{T}}\f$,
 * \f$f_{\text{opt}} = 3\f$
 * (domain scaled to \f$[0, 1]^2\f$)
 */
class GoldsteinPrice : public Test
{
public:
    /**
     * Constructor.
     */
    GoldsteinPrice() : Test(2)
    {
    }
    
    /**
     * Evaluates the test function.
     * 
     * @param x     point \f$\vec{x} \in [0, 1]^2\f$
     * @return      \f$f(\vec{x})\f$
     */
    double evalUndisplaced(const std::vector<double> &x)
    {
        const double x1 = 4.0 * x[0] - 2.0;
        const double x2 = 4.0 * x[1] - 2.0;
        
        return (1.0 + (x1+x2+1.0)*(x1+x2+1.0) *
                      (19.0 - 14.0*x1 + 3.0*x1*x1 - 14.0*x2 + 6.0*x1*x2 + 3.0*x2*x2)) *
               (30.0 + (2.0*x1-3.0*x2)*(2.0*x1-3.0*x2) *
                       (18.0 - 32.0*x1 + 12.0*x1*x1 + 48.0*x2 - 36.0*x1*x2 + 27.0*x2*x2));
    }
    
    /**
     * Returns minimal point and function value of the test function.
     * 
     * @param[out] x    minimal point \f$\vec{x}_{\text{opt}} \in [0, 1]^2\f$
     * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
     */
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x.clear();
        x.push_back(0.5);
        x.push_back(0.25);
        return 3.0;
    }
    
    /**
     * @return clone of the object
     */
    virtual tools::SmartPointer<Objective> clone()
    {
        return tools::SmartPointer<Objective>(new GoldsteinPrice(*this));
    }
};

}
}
}
}

#endif
