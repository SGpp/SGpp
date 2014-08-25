/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_OPERATION_OPERATIONEVALHESSIAN_HPP
#define SGPP_OPT_OPERATION_OPERATIONEVALHESSIAN_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

namespace sg
{
namespace opt
{

/**
 * Abstract operation for evaluating a linear combination of basis functions, its gradient
 * and its Hessian.
 */
class OperationEvalHessian
{
public:
    /**
     * Constructor.
     */
    OperationEvalHessian()
    {
    }
    
    /**
     * Virtual destructor.
     */
    virtual ~OperationEvalHessian()
    {
    }
    
    /**
     * Pure virtual method for evaluating a linear combination of basis functions, its gradient
     * and its Hessian.
     * 
     * @param       alpha       coefficient vector
     * @param       point       evaluation point
     * @param[out]  gradient    gradient vector of the linear combination
     * @param[out]  hessian     Hessian matrix of the linear combination
     * @return                  value of the linear combination
     */
    virtual double evalHessian(base::DataVector &alpha, const std::vector<double> &point,
                               base::DataVector &gradient, base::DataMatrix &hessian) = 0;
    
    /**
     * Convenience function for using base::DataVector as points.
     * 
     * @param       alpha       coefficient vector
     * @param       point       evaluation point
     * @param[out]  gradient    gradient vector of the linear combination
     * @param[out]  hessian     Hessian matrix of the linear combination
     * @return                  value of the linear combination
     */
    virtual double evalHessian(base::DataVector &alpha, base::DataVector &point,
                               base::DataVector &gradient, base::DataMatrix &hessian)
    {
        const std::vector<double> p(point.getPointer(), point.getPointer() + point.getSize());
        return evalHessian(alpha, p, gradient, hessian);
    }
};

}
}

#endif
