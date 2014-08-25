/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_FUNCTION_INTERPOLANT_HPP
#define SGPP_OPT_FUNCTION_INTERPOLANT_HPP

#include "opt/function/Objective.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/operation/OperationEval.hpp"
#include "base/grid/Grid.hpp"
#include "opt/operation/OpFactory.hpp"

#include <vector>
#include <cstring>

namespace sg
{
namespace opt
{
namespace function
{

/**
 * Sparse grid interpolant as an objective function.
 * 
 * More generally, the function can be any linear combination
 * \f$f\colon [0, 1]^d \to \mathbb{R}\f$,
 * \f$f(\vec{x}) = \sum_{k=1}^N \alpha_k \varphi_k(\vec{x})\f$ of the basis functions
 * \f$\varphi_k = \varphi_{\vec{\ell}_k,\vec{i}_k}\f$ of a sparse grid with grid points
 * \f$\vec{x}_k = \vec{x}_{\vec{\ell}_k,\vec{i}_k}\f$.
 * But most often, the function (e.g. its coefficients) is constructed as an interpolant
 * at the grid points for some function values.
 */
class Interpolant : public Objective
{
public:
    /**
     * Constructor.
     * Do not destruct the grid before the Interpolant object!
     * 
     * @param d     dimension of the domain
     * @param grid  sparse grid
     * @param alpha coefficient vector
     */
    Interpolant(size_t d, base::Grid &grid, base::DataVector &alpha) :
        Objective(d),
        grid(grid),
        op_eval(createOperationEval(grid)),
        alpha(alpha)
    {
    }
    
    /**
     * Evaluation of the function.
     * 
     * @param x     point \f$\vec{x} \in \mathbb{R}^d\f$
     * @return      \f$f(\vec{x})\f$
     */
    inline double eval(const std::vector<double> &x)
    {
        // copy x, necessary due to non-existing const correctness in sg::base
        std::vector<double> y = x;
        return op_eval->eval(alpha, y);
    }
    
    /**
     * @return clone of the object
     */
    virtual tools::SmartPointer<Objective> clone()
    {
        return tools::SmartPointer<Objective>(new Interpolant(d, grid, alpha));
    }
    
protected:
    /// sparse grid
    base::Grid &grid;
    /// pointer to evaluation operation
    tools::SmartPointer<base::OperationEval> op_eval;
    /// coefficient vector
    base::DataVector alpha;
};

}
}
}

#endif
