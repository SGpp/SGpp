/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_TEST_SPHERE_HPP
#define SGPP_OPT_FUNCTION_TEST_SPHERE_HPP

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
 * Sphere test function.
 * 
 * Definition:
 * \f$f(\vec{x}) := \lVert \vec{x} \rVert_2^2\f$,
 * \f$\vec{x} \in [-1, 9]^d\f$,
 * \f$\vec{x}_{\text{opt}} = \vec{0}\f$,
 * \f$f_{\text{opt}} = 0\f$
 * (domain scaled to \f$[0, 1]^d\f$)
 */
class Sphere : public Test
{
public:
    /**
     * Constructor.
     * 
     * @param d     dimension of the domain
     */
    Sphere(size_t d) : Test(d)
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
        double result = 0.0;
        
        for (size_t t = 0; t < d; t++)
        {
            const double xt = 10.0 * x[t] - 1.0;
            result += xt*xt;
        }
        
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
        x = std::vector<double>(d, 0.1);
        return 0.0;
    }
    
    /**
     * @return clone of the object
     */
    virtual tools::SmartPointer<Objective> clone()
    {
        return tools::SmartPointer<Objective>(new Sphere(*this));
    }
};

}
}
}
}

#endif
