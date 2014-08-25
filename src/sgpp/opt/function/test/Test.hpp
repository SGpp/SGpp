/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_TEST_TEST_HPP
#define SGPP_OPT_FUNCTION_TEST_TEST_HPP

#include <vector>
#include <cstddef>

#include "opt/function/Objective.hpp"

namespace sg
{
namespace opt
{
namespace function
{
namespace test
{

/**
 * Base class for analytical objective function examples ("test functions").
 * The only difference to Objective is the possibility for pseudorandom displacements of the
 * function and the specification of the (or an) minimal point.
 * 
 * Taking the average of results of multiple runs with different displacements makes
 * results more robust and significant.
 * The displaced function is \f$\vec{x} \mapsto f(\vec{x} + \vec{d})\f$ for a vector \f$\vec{d}\f$
 * ("displacement") with each component \f$d_t\f$ distributed normally with mean 0 and
 * a specific standard deviation.
 */
class Test : public Objective
{
public:
    /// default standard deviation
    static const double DEFAULT_STANDARD_DEVIATION;
    
    /**
     * Constructor.
     * The displacement is set to all zeros, so to displace the function call
     * generateDisplacement() afterwards.
     * 
     * @param d     dimension of the domain
     */
    Test(size_t d);
    
    /**
     * Virtual destructor.
     */
    virtual ~Test();
    
    /**
     * Evaluate displaced function.
     * 
     * @param x     point \f$\vec{x} \in \mathbb{R}^d\f$
     * @return      \f$f(\vec{x} + \vec{d})\f$ with displacement \f$\vec{d}\f$
     */
    double eval(const std::vector<double> &x);
    
    /**
     * Pure virtual method for evaluating the undisplaced function.
     * 
     * @param x     point \f$\vec{x} \in \mathbb{R}^d\f$
     * @return      \f$f(\vec{x})\f$
     */
    virtual double evalUndisplaced(const std::vector<double> &x) = 0;
    
    /**
     * Returns the minimal point of the displaced function.
     * 
     * @param[out] x    reverse displaced minimal point \f$\vec{x}_{\text{opt}} - \vec{d}\f$
     * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
     */
    double getOptimalPoint(std::vector<double> &x);
    
    /**
     * Pure virtual method returning the minimal point (or one of the minimal points, if
     * there are multiple of them) of the test function.
     * 
     * @param[out] x    minimal point \f$\vec{x}_{\text{opt}}\f$
     * @return          minimal function value \f$f_{\text{opt}} = f(\vec{x}_{\text{opt}})\f$
     */
    virtual double getOptimalPointUndisplaced(std::vector<double> &x) = 0;
    
    /**
     * Generate normally distributed pseudorandom displacement with default standard deviation.
     */
    void generateDisplacement();
    
    /**
     * Generate normally distributed pseudorandom displacement.
     * This function can be overridden if the minimal points of the test function lie
     * near or on the boundary or if there is the chance of a pole getting in the domain
     * by displacing.
     * 
     * @param std_dev   standard deviation of the displacement coordinates
     */
    virtual void generateDisplacement(double std_dev);
    
    /**
     * Add the displacement to a vector.
     * 
     * @param[in,out] x     vector to be displaced
     */
    void displaceVector(std::vector<double> &x) const;
    
    /**
     * Subtract the displacement from a vector.
     * 
     * @param[in,out] x     vector to be reverse displaced
     */
    void reverseDisplaceVector(std::vector<double> &x) const;
    
    /**
     * @return standard deviation of the displacement
     */
    double getStandardDeviation() const;
    
    /**
     * @param[out] displacement     currently used displacement
     */
    void getDisplacement(std::vector<double> &displacement) const;
    
protected:
    /// standard deviation
    double std_dev;
    /// vector displacement
    std::vector<double> displacement;
};

}
}
}
}

#endif
