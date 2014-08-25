/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_TEST_EASOM_HPP
#define SGPP_OPT_FUNCTION_TEST_EASOM_HPP

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
 * Easom test function.
 * 
 * Definition:
 * \f$f(\vec{x}) := -\cos x_1 \cos x_2 \exp(-(x_1 - \pi)^2 - (x_2 - \pi)^2)\f$,
 * \f$\vec{x} \in [-100, 100]^2\f$,
 * \f$\vec{x}_{\text{opt}} = (\pi, \pi)^{\mathrm{T}}\f$,
 * \f$f_{\text{opt}} = -1\f$
 * (domain scaled to \f$[0, 1]^2\f$)
 */
class Easom : public Test
{
public:
    /**
     * Constructor.
     */
    Easom() : Test(2)
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
        const double x1 = 200.0 * x[0] - 100.0;
        const double x2 = 200.0 * x[1] - 100.0;
        
        return -std::cos(x1) * std::cos(x2) *
                std::exp(-((x1-M_PI)*(x1-M_PI) + (x2-M_PI)*(x2-M_PI)));
    }
    
    /**
     * Returns minimal point and function value of the test function.
     * 
     * @param[out] x    minimal point \f$\vec{x}_{\text{opt}} \in [0, 1]^2\f$
     * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
     */
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x = std::vector<double>(2, 0.51570796326794896619231);
        return -1.0;
    }
    
    /**
     * @return clone of the object
     */
    virtual tools::SmartPointer<Objective> clone()
    {
        return tools::SmartPointer<Objective>(new Easom(*this));
    }
};

}
}
}
}

#endif
