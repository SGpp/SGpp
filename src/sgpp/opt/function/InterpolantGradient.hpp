/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_INTERPOLANTGRADIENT_HPP
#define SGPP_OPT_FUNCTION_INTERPOLANTGRADIENT_HPP

#include "opt/function/ObjectiveGradient.hpp"
#include "base/datatypes/DataVector.hpp"
#include "opt/operation/OperationEvalGradient.hpp"
#include "opt/operation/OpFactory.hpp"

#include <vector>

namespace sg
{
namespace opt
{
namespace function
{

/**
 * Sparse grid interpolant gradient as an objective gradient.
 * 
 * @see Interpolant
 */
class InterpolantGradient : public ObjectiveGradient
{
public:
    /**
     * Constructor.
     * Do not destruct the grid before the InterpolantGradient object!
     * 
     * @param d     dimension of the domain
     * @param grid  sparse grid
     * @param alpha coefficient vector
     */
    InterpolantGradient(size_t d, base::Grid &grid, base::DataVector &alpha) :
            ObjectiveGradient(d), grid(grid),
            op_eval_gradient(createOperationEvalGradient(grid)),
            alpha(alpha)
    {
    }
    
    /**
     * Evaluation of the function and its gradient.
     * 
     * @param      x            point \f$\vec{x} \in \mathbb{R}^d\f$
     * @param[out] gradient     gradient \f$\nabla f(\vec{x}) \in \mathbb{R}^d\f$
     * @return                  \f$f(\vec{x})\f$
     */
    inline double evalGradient(const std::vector<double> &x, base::DataVector &gradient)
    {
        // copy x, necessary due to non-existing const correctness in sg::base
        std::vector<double> y = x;
        return op_eval_gradient->evalGradient(alpha, y, gradient);
    }
    
    /**
     * @return clone of the object
     */
    virtual tools::SmartPointer<ObjectiveGradient> clone()
    {
        return tools::SmartPointer<ObjectiveGradient>(new InterpolantGradient(d, grid, alpha));
    }
    
protected:
    /// sparse grid
    base::Grid &grid;
    /// pointer to evaluation operation
    tools::SmartPointer<OperationEvalGradient> op_eval_gradient;
    /// coefficient vector
    base::DataVector alpha;
};

}
}
}

#endif
