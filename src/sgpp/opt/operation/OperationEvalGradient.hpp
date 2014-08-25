/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_OPERATION_OPERATIONEVALGRADIENT_HPP
#define SGPP_OPT_OPERATION_OPERATIONEVALGRADIENT_HPP

#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace opt
{

/**
 * Abstract operation for evaluating a linear combination of basis functions and its gradient.
 */
class OperationEvalGradient
{
public:
    /**
     * Constructor.
     */
    OperationEvalGradient()
    {
    }
    
    /**
     * Virtual destructor.
     */
    virtual ~OperationEvalGradient()
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
