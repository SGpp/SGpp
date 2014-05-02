/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALBSPLINEBOUNDARY_HPP
#define OPERATIONEVALBSPLINEBOUNDARY_HPP

#include "base/operation/OperationEval.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/bspline/boundary/BsplineBoundaryBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace base
{

class OperationEvalBsplineBoundary : public OperationEval
{
public:
    OperationEvalBsplineBoundary(GridStorage *storage, size_t degree) :
            storage(storage), base(degree) {}
    
    virtual ~OperationEvalBsplineBoundary() {}
    
    virtual double eval(DataVector &alpha, std::vector<double> &point);
    
protected:
    GridStorage *storage;
    SBsplineBoundaryBase base;
};

}
}

#endif /* OPERATIONEVALBSPLINEBOUNDARY_HPP */
