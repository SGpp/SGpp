/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALGRADIENTMODWAVELET_HPP
#define OPERATIONEVALGRADIENTMODWAVELET_HPP

#include "base/operation/OperationEvalGradient.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/modwavelet/ModifiedWaveletBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace base
{

class OperationEvalGradientModWavelet : public OperationEvalGradient
{
public:
    OperationEvalGradientModWavelet(GridStorage *storage) : storage(storage) {}
    virtual ~OperationEvalGradientModWavelet() {}
    
    virtual double evalGradient(DataVector &alpha, const std::vector<double> &point,
                                DataVector &gradient);
    
protected:
    GridStorage* storage;
    SModWaveletBase base;
};

}
}

#endif /* OPERATIONEVALGRADIENTMODWAVELET_HPP */
