/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALGRADIENTBSPLINEBOUNDARY_HPP
#define OPERATIONEVALGRADIENTBSPLINEBOUNDARY_HPP

#include "base/operation/OperationEvalGradient.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/bspline/boundary/BsplineBoundaryBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace base
{

class OperationEvalGradientBsplineBoundary : public OperationEvalGradient
{
public:
    OperationEvalGradientBsplineBoundary(GridStorage *storage, size_t degree) :
            storage(storage), base(degree) {}
    
    virtual ~OperationEvalGradientBsplineBoundary() {}
    
    virtual double evalGradient(DataVector &alpha, const std::vector<double> &point,
                                DataVector &gradient);
    
protected:
    GridStorage *storage;
    SBsplineBoundaryBase base;
};

}
}

#endif /* OPERATIONEVALGRADIENTBSPLINEBOUNDARY_HPP */
