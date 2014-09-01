/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_BASE_OPERATION_OPERATIONNAIVEEVALGRADIENT_HPP
#define SGPP_BASE_OPERATION_OPERATIONNAIVEEVALGRADIENT_HPP

#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace base
{

/**
 * Abstract operation for evaluating a linear combination of basis functions and its gradient.
 * The "naive" is indicating that classes implementing this operation should use a "naive"
 * approach, e.g. by evaluating all basis functions by brute force.
 */
class OperationNaiveEvalGradient
{
public:
    /**
     * Constructor.
     */
    OperationNaiveEvalGradient()
    {
    }
    
    /**
     * Virtual destructor.
     */
    virtual ~OperationNaiveEvalGradient()
    {
    }
    
    /**
     * Pure virtual method for evaluating a linear combination of basis functions and its gradient.
     * 
     * @param       alpha       coefficient vector
     * @param       point       evaluation point
     * @param[out]  gradient    gradient vector of the linear combination
     * @return                  value of the linear combination
     */
    virtual double evalGradient(base::DataVector &alpha, const std::vector<double> &point,
                                base::DataVector &gradient) = 0;
    
    /**
     * Convenience function for using base::DataVector as points.
     * 
     * @param       alpha       coefficient vector
     * @param       point       evaluation point
     * @param[out]  gradient    gradient vector of the linear combination
     * @return                  value of the linear combination
     */
    virtual double evalGradient(base::DataVector &alpha, base::DataVector &point,
                                base::DataVector &gradient)
    {
        const std::vector<double> p(point.getPointer(), point.getPointer() + point.getSize());
        return evalGradient(alpha, p, gradient);
    }
};

}
}

#endif
