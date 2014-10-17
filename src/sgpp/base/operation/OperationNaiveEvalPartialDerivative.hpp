/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_BASE_OPERATION_OPERATIONNAIVEEVALPARTIALDERIVATIVE_HPP
#define SGPP_BASE_OPERATION_OPERATIONNAIVEEVALPARTIALDERIVATIVE_HPP

#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace base
{

/**
 * Abstract operation for evaluating a partial derivative of a linear combination of
 * basis functions.
 * The "naive" is indicating that classes implementing this operation should use a "naive"
 * approach, e.g. by evaluating all basis functions by brute force.
 */
class OperationNaiveEvalPartialDerivative
{
public:
    /**
     * Constructor.
     */
    OperationNaiveEvalPartialDerivative()
    {
    }
    
    /**
     * Virtual destructor.
     */
    virtual ~OperationNaiveEvalPartialDerivative()
    {
    }
    
    /**
     * Pure virtual method for evaluating a partial derivative of a linear combination of
     * basis functions.
     * 
     * @param alpha     coefficient vector
     * @param point     evaluation point
     * @param deriv_dim dimension in which the partial derivative should be taken (0, ..., d-1)
     * @return          value of the partial derivative of the linear combination
     */
    virtual double evalPartialDerivative(base::DataVector& alpha,
                                         const std::vector<double>& point,
                                         size_t deriv_dim) = 0;
    
    /**
     * Convenience function for using base::DataVector as points.
     * 
     * @param alpha     coefficient vector
     * @param point     evaluation point
     * @param deriv_dim dimension in which the partial derivative should be taken (0, ..., d-1)
     * @return          value of the partial derivative of the linear combination
     */
    virtual double evalPartialDerivative(base::DataVector& alpha,
                                         base::DataVector& point,
                                         size_t deriv_dim)
    {
        const std::vector<double> p(point.getPointer(), point.getPointer() + point.getSize());
        return evalPartialDerivative(alpha, p, deriv_dim);
    }
};

}
}

#endif
