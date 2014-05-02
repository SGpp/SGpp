/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALWAVELET_HPP
#define OPERATIONEVALWAVELET_HPP

#include "base/operation/OperationEval.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/wavelet/noboundary/WaveletBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace base
{

class OperationEvalWavelet : public OperationEval
{
public:
    OperationEvalWavelet(GridStorage* storage) : storage(storage) {}
    virtual ~OperationEvalWavelet() {}
    
    virtual double eval(DataVector &alpha, std::vector<double> &point);
    
protected:
    GridStorage *storage;
    SWaveletBase base;
};

}
}

#endif /* OPERATIONEVALWAVELET_HPP */
