/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALGRADIENTWAVELET_HPP
#define OPERATIONEVALGRADIENTWAVELET_HPP

#include "base/operation/OperationEvalGradient.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/wavelet/noboundary/WaveletBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace base
{

class OperationEvalGradientWavelet : public OperationEvalGradient
{
public:
    OperationEvalGradientWavelet(GridStorage *storage) : storage(storage) {}
    virtual ~OperationEvalGradientWavelet() {}
    
    virtual double evalGradient(DataVector &alpha, const std::vector<double> &point,
                                DataVector &gradient);
    
protected:
    GridStorage* storage;
    SWaveletBase base;
};

}
}

#endif /* OPERATIONEVALGRADIENTWAVELET_HPP */
