/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_TEST_HARTMAN3_HPP
#define SGPP_OPT_FUNCTION_TEST_HARTMAN3_HPP

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
 * Hartman3 test function.
 * 
 * Definition:
 * \f$f(\vec{x}) := -\sum_{i=1}^4 a_i
 *                  \exp\!\left(-\sum_{t=1}^3 b_{i,t} (x_t - c_{i,t})^2\right)\f$,
 * \f$\vec{a} = \begin{pmatrix}1\\1.2\\3\\3.2\end{pmatrix}\f$,
 * \f$B := 
 *      \begin{pmatrix}
 *          3 & 10 & 30\\
 *          0.1 & 10 & 35\\
 *          3 & 10 & 30\\
 *          0.1 & 10 & 35
 *      \end{pmatrix}\f$,
 * \f$C :=
 *      \begin{pmatrix}
 *          0.3689 & 0.1170 & 0.2673\\
 *          0.4699 & 0.4387 & 0.7470\\
 *          0.1091 & 0.8732 & 0.5547\\
 *          0.0382 & 0.5743 & 0.8828
 *      \end{pmatrix}\f$,
 * \f$\vec{x} \in [0, 1]^3\f$,
 * \f$\vec{x}_{\text{opt}} = (0.114614, 0.555649, 0.852547)^{\mathrm{T}}\f$,
 * \f$f_{\text{opt}} = -3.862785\f$
 */
class Hartman3 : public Test
{
public:
    /**
     * Constructor.
     */
    Hartman3() : Test(3)
    {
    }
    
    /**
     * Evaluates the test function.
     * 
     * @param x     point \f$\vec{x} \in [0, 1]^3\f$
     * @return      \f$f(\vec{x})\f$
     */
    double evalUndisplaced(const std::vector<double> &x)
    {
        return -1.0 * exp(-3.0*(x[0]-0.3689)*(x[0]-0.3689) - 10.0*(x[1]-0.1170)*(x[1]-0.1170) -
                          30.0*(x[2]-0.2673)*(x[2]-0.2673)) -
                1.2 * exp(-0.1*(x[0]-0.4699)*(x[0]-0.4699) - 10.0*(x[1]-0.4387)*(x[1]-0.4387) -
                          35.0*(x[2]-0.7470)*(x[2]-0.7470)) -
                3.0 * exp(-3.0*(x[0]-0.1091)*(x[0]-0.1091) - 10.0*(x[1]-0.8732)*(x[1]-0.8732) -
                          30.0*(x[2]-0.5547)*(x[2]-0.5547)) -
                3.2 * exp(-0.1*(x[0]-0.0382)*(x[0]-0.0382) - 10.0*(x[1]-0.5743)*(x[1]-0.5743) -
                          35.0*(x[2]-0.8828)*(x[2]-0.8828));
    }
    
    /**
     * Returns minimal point and function value of the test function.
     * 
     * @param[out] x    minimal point \f$\vec{x}_{\text{opt}} \in [0, 1]^3\f$
     * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
     */
    double getOptimalPointUndisplaced(std::vector<double> &x)
    {
        x.clear();
        x.push_back(0.114614);
        x.push_back(0.555649);
        x.push_back(0.852547);
        return evalUndisplaced(x);
    }
    
    /**
     * @return clone of the object
     */
    virtual tools::SmartPointer<Objective> clone()
    {
        return tools::SmartPointer<Objective>(new Hartman3(*this));
    }
};

}
}
}
}

#endif
