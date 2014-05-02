/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALHESSIANWAVELET_HPP
#define OPERATIONEVALHESSIANWAVELET_HPP

#include "base/operation/OperationEvalHessian.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/wavelet/noboundary/WaveletBasis.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

namespace sg
{
namespace base
{

class OperationEvalHessianWavelet : public OperationEvalHessian
{
public:
    OperationEvalHessianWavelet(GridStorage *storage) : storage(storage) {}
    virtual ~OperationEvalHessianWavelet() {}
    
    virtual double evalHessian(DataVector &alpha, const std::vector<double> &point,
                               DataVector &gradient, DataMatrix &hessian);
    
protected:
    GridStorage* storage;
    SWaveletBase base;
};

}
}

#endif /* OPERATIONEVALHESSIANWAVELET_HPP */
