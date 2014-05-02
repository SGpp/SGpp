/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALHESSIANMODWAVELET_HPP
#define OPERATIONEVALHESSIANMODWAVELET_HPP

#include "base/operation/OperationEvalHessian.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/modwavelet/ModifiedWaveletBasis.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

namespace sg
{
namespace base
{

class OperationEvalHessianModWavelet : public OperationEvalHessian
{
public:
    OperationEvalHessianModWavelet(GridStorage *storage) : storage(storage) {}
    virtual ~OperationEvalHessianModWavelet() {}
    
    virtual double evalHessian(DataVector &alpha, const std::vector<double> &point,
                               DataVector &gradient, DataMatrix &hessian);
    
protected:
    GridStorage* storage;
    SModWaveletBase base;
};

}
}

#endif /* OPERATIONEVALHESSIANMODWAVELET_HPP */
