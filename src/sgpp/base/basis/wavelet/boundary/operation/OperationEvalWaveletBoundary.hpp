/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALWAVELETBOUNDARY_HPP
#define OPERATIONEVALWAVELETBOUNDARY_HPP

#include "base/operation/OperationEval.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/wavelet/boundary/WaveletBoundaryBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace base
{

class OperationEvalWaveletBoundary : public OperationEval
{
public:
    OperationEvalWaveletBoundary(GridStorage* storage) : storage(storage) {}
    virtual ~OperationEvalWaveletBoundary() {}
    
    virtual double eval(DataVector &alpha, std::vector<double> &point);
    
protected:
    GridStorage *storage;
    SWaveletBoundaryBase base;
};

}
}

#endif /* OPERATIONEVALWAVELETBOUNDARY_HPP */
