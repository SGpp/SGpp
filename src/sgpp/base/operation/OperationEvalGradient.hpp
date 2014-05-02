/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALGRADIENT_HPP
#define OPERATIONEVALGRADIENT_HPP

#include "base/datatypes/DataVector.hpp"

#ifdef _WIN32
#pragma warning(disable: 4267)
#endif

namespace sg
{
namespace base
{

class OperationEvalGradient
{
public:
    OperationEvalGradient() {}
    virtual ~OperationEvalGradient() {}
    
    virtual double evalGradient(DataVector &alpha, const std::vector<double> &point,
                                DataVector &gradient) = 0;
    
    double evalGradient(DataVector &alpha, DataVector &point, DataVector &gradient)
    {
        std::vector<double> p(point.getPointer(), point.getPointer() + point.getSize());
        return evalGradient(alpha, p, gradient);
    }
};

}
}

#endif /* OPERATIONEVALGRADIENT_HPP */
