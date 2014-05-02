/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALBSPLINECLENSHAWCURTIS_HPP
#define OPERATIONEVALBSPLINECLENSHAWCURTIS_HPP

#include "base/operation/OperationEval.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/bspline/clenshawcurtis/BsplineClenshawCurtisBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace base
{

class OperationEvalBsplineClenshawCurtis : public OperationEval
{
public:
    OperationEvalBsplineClenshawCurtis(GridStorage *storage, size_t degree) :
            OperationEvalBsplineClenshawCurtis(storage, degree, NULL) {}
    OperationEvalBsplineClenshawCurtis(GridStorage *storage, size_t degree,
                                       const CosineTable *cosine_table) :
            storage(storage),
            base(degree, cosine_table) {}
    
    virtual ~OperationEvalBsplineClenshawCurtis() {}
    
    virtual double eval(DataVector &alpha, std::vector<double> &point);
    
protected:
    GridStorage *storage;
    SBsplineClenshawCurtisBase base;
};

}
}

#endif /* OPERATIONEVALBSPLINECLENSHAWCURTIS_HPP */
