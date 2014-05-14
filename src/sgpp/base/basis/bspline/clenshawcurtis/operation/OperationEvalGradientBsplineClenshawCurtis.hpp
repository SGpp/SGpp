/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALGRADIENTBSPLINECLENSHAWCURTIS_HPP
#define OPERATIONEVALGRADIENTBSPLINECLENSHAWCURTIS_HPP

#include "base/operation/OperationEvalGradient.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/bspline/clenshawcurtis/BsplineClenshawCurtisBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace base
{

class OperationEvalGradientBsplineClenshawCurtis : public OperationEvalGradient
{
public:
    OperationEvalGradientBsplineClenshawCurtis(GridStorage *storage, size_t degree,
                                               const CosineTable *cosine_table = nullptr) :
            storage(storage),
            base(degree, cosine_table) {}
    
    virtual ~OperationEvalGradientBsplineClenshawCurtis() {}
    
    virtual double evalGradient(DataVector &alpha, const std::vector<double> &point,
                                DataVector &gradient);
    
protected:
    GridStorage *storage;
    SBsplineClenshawCurtisBase base;
};

}
}

#endif /* OPERATIONEVALGRADIENTBSPLINECLENSHAWCURTIS_HPP */
