/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALHESSIAN_HPP
#define OPERATIONEVALHESSIAN_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

#ifdef _WIN32
#pragma warning(disable: 4267)
#endif

namespace sg
{
namespace base
{

class OperationEvalHessian
{
public:
    OperationEvalHessian() {}
    virtual ~OperationEvalHessian() {}
    
    virtual double evalHessian(DataVector &alpha, const std::vector<double> &point,
                             DataVector &gradient, DataMatrix &hessian) = 0;
    
    double evalHessian(DataVector &alpha, DataVector &point,
                       DataVector &gradient, DataMatrix &hessian)
    {
        std::vector<double> p(point.getPointer(), point.getPointer() + point.getSize());
        return evalHessian(alpha, p, gradient, hessian);
    }
};

}
}

#endif /* OPERATIONEVALHESSIAN_HPP */
