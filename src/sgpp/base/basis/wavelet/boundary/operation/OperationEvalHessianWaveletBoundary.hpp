/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALHESSIANBOUNDARYWAVELET_HPP
#define OPERATIONEVALHESSIANBOUNDARYWAVELET_HPP

#include "base/operation/OperationEvalHessian.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/wavelet/boundary/WaveletBoundaryBasis.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

namespace sg
{
namespace base
{

class OperationEvalHessianWaveletBoundary : public OperationEvalHessian
{
public:
    OperationEvalHessianWaveletBoundary(GridStorage *storage) : storage(storage) {}
    virtual ~OperationEvalHessianWaveletBoundary() {}
    
    virtual double evalHessian(DataVector &alpha, const std::vector<double> &point,
                               DataVector &gradient, DataMatrix &hessian);
    
protected:
    GridStorage* storage;
    SWaveletBoundaryBase base;
};

}
}

#endif /* OPERATIONEVALHESSIANBOUNDARYWAVELET_HPP */
