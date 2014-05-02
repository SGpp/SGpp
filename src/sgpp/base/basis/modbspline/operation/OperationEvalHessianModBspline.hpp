/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALHESSIANMODBSPLINE_HPP
#define OPERATIONEVALHESSIANMODBSPLINE_HPP

#include "base/operation/OperationEvalHessian.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/modbspline/ModifiedBsplineBasis.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

namespace sg
{
namespace base
{

class OperationEvalHessianModBspline : public OperationEvalHessian
{
public:
    OperationEvalHessianModBspline(GridStorage *storage, size_t degree) :
            storage(storage), base(degree) {}
    
    virtual ~OperationEvalHessianModBspline() {}
    
    virtual double evalHessian(DataVector &alpha, const std::vector<double> &point,
                               DataVector &gradient, DataMatrix &hessian);
    
protected:
    GridStorage *storage;
    SModBsplineBase base;
};

}
}

#endif /* OPERATIONEVALHESSIANMODBSPLINE_HPP */
