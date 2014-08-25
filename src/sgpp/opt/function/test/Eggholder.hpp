/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_TEST_EGGHOLDER_HPP
#define SGPP_OPT_FUNCTION_TEST_EGGHOLDER_HPP

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
 * Eggholder test function.
 * 
 * Definition:
 * \f$f(\vec{x}) := -(x_2 + 47) \sin\!\left(\sqrt{\left|x_1/2 + x_2 + 47\right|}\right)
 *                - x_1 \sin\!\left(\sqrt{\left|x_1 - x_2 - 47\right|}\right)\f$,
 * \f$\vec{x} \in [-512, 512]^2\f$,
 * \f$\vec{x}_{\text{opt}} = (512, 404.2319)^{\mathrm{T}}\f$,
 * \f$f_{\text{opt}} = -959.6407\f$
 * (domain scaled to \f$[0, 1]^2\f$)
 * 
 * The displacement is restricted because the minimal point lies on the boundary of \f$[0, 1]^2\f$.
 */
class Eggholder : public Test
{
public:
    /**
     * Constructor.
     */
    Eggholder() : Test(2)
    {
    }
    
    /**
     * Generate normally distributed pseudorandom displacement
     * with the restriction of \f$d_1 = 0\f$.
     * 
     * @param std_dev   standard deviation of the displacement coordinates
     */
    void generateDisplacement(double std_dev)
    {
        Test::generateDisplacement(std_dev);
        displacement[0] = 0.0;
    }
    
    /**
     * Evaluates the test function.
     * 
     * @param x     point \f$\vec{x} \in [0, 1]^2\f$
     * @return      \f$f(\vec{x})\f$
     */
    double evalUndisplaced(const std::vector<double> &x)
    {
        const double x1 = 1024.0 * x[0] - 512.0;
        const double x2 = 1024.0 * x[1] - 512.0;
        
        return -(x2 + 47.0) * std::sin(std::sqrt(std::abs(x1/2.0 + x2 + 47.0))) -
                x1 * std::sin(std::sqrt(std::abs(x1 - (x2 + 47.0))));
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
        x.push_back(1.0);
        x.push_back(0.8947577);
        return evalUndisplaced(x);
    }
    
    /**
     * @return clone of the object
     */
    virtual tools::SmartPointer<Objective> clone()
    {
        return tools::SmartPointer<Objective>(new Eggholder(*this));
    }
};

}
}
}
}

#endif
