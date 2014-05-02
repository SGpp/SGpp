/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALBSPLINE_HPP
#define OPERATIONEVALBSPLINE_HPP

#include "base/operation/OperationEval.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/bspline/noboundary/BsplineBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace base
{

class OperationEvalBspline : public OperationEval
{
public:
    OperationEvalBspline(GridStorage *storage, size_t degree) :
            storage(storage), base(degree) {}
    
    virtual ~OperationEvalBspline() {}
    
    virtual double eval(DataVector &alpha, std::vector<double> &point);
    
protected:
    GridStorage *storage;
    SBsplineBase base;
};

}
}

#endif /* OPERATIONEVALBSPLINE_HPP */
