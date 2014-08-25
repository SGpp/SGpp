/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_TEST_MLADINEO_HPP
#define SGPP_OPT_FUNCTION_TEST_MLADINEO_HPP

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
 * Mladineo test function.
 * 
 * Definition:
 * \f$f(\vec{x}) := 1 + \frac{1}{2} \lVert \vec{x} \rVert_2^2
 *                    - \cos(10 \ln(2x_1)) \cos(10 \ln(3x_2))\f$,
 * \f$\vec{x} \in [0.01, 1]^2\f$,
 * \f$\vec{x}_{\text{opt}} = (0.0115, 0.0144)^{\mathrm{T}}\}\f$,
 * \f$f_{\text{opt}} = 0.000170\f$
 * (domain scaled to \f$[0, 1]^2\f$)
 * 
 * The displacement is restricted because the minimal points lie near the origin.
 */
class Mladineo : public Test
{
public:
    /**
     * Constructor.
     */
    Mladineo() : Test(2)
    {
    }
    
    /**
     * Generate normally distributed pseudorandom displacement
     * with the restriction of \f$\vec{d} \in [-0.01, 0] \times [-0.01, 0]\f$.
     * 
     * @param std_dev   standard deviation of the displacement coordinates
     */
    void generateDisplacement(double std_dev)
    {
        do
        {
            Test::generateDisplacement(std_dev);
        } while ((displacement[0] > 0) || (displacement[0] < -0.01) ||
                 (displacement[1] > 0) || (displacement[1] < -0.01));
    }
    
    /**
     * Evaluates the test function.
     * 
     * @param x     point \f$\vec{x} \in [0, 1]^2\f$
     * @return      \f$f(\vec{x})\f$
     */
    double evalUndisplaced(const std::vector<double> &x)
    {
        const double x1 = 0.99 * x[0] + 0.01;
        const double x2 = 0.99 * x[1] + 0.01;
        
        return 1.0 + (x1*x1 + x2*x2) / 2.0 - cos(10.0*log(2.0*x1))*cos(10.0*log(3.0*x2));
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
        x.push_back(0.001542);
        x.push_back(0.004449);
        return evalUndisplaced(x);
    }
    
    /**
     * @return clone of the object
     */
    virtual tools::SmartPointer<Objective> clone()
    {
        return tools::SmartPointer<Objective>(new Mladineo(*this));
    }
};

}
}
}
}

#endif
