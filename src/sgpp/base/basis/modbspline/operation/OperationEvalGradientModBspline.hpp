/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALGRADIENTMODBSPLINE_HPP
#define OPERATIONEVALGRADIENTMODBSPLINE_HPP

#include "base/operation/OperationEvalGradient.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/modbspline/ModifiedBsplineBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace base
{

class OperationEvalGradientModBspline : public OperationEvalGradient
{
public:
    OperationEvalGradientModBspline(GridStorage *storage, size_t degree) :
            storage(storage), base(degree) {}
    
    virtual ~OperationEvalGradientModBspline() {}
    
    virtual double evalGradient(DataVector &alpha, const std::vector<double> &point,
                                DataVector &gradient);
    
protected:
    GridStorage *storage;
    SModBsplineBase base;
};

}
}

#endif /* OPERATIONEVALGRADIENTMODBSPLINE_HPP */
