/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALHESSIANBSPLINECLENSHAWCURTIS_HPP
#define OPERATIONEVALHESSIANBSPLINECLENSHAWCURTIS_HPP

#include "base/operation/OperationEvalHessian.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/bspline/clenshawcurtis/BsplineClenshawCurtisBasis.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

namespace sg
{
namespace base
{

class OperationEvalHessianBsplineClenshawCurtis : public OperationEvalHessian
{
public:
    OperationEvalHessianBsplineClenshawCurtis(GridStorage *storage, size_t degree) :
            OperationEvalHessianBsplineClenshawCurtis(storage, degree, NULL) {}
    OperationEvalHessianBsplineClenshawCurtis(GridStorage *storage, size_t degree,
                                              const CosineTable *cosine_table) :
            storage(storage),
            base(degree, cosine_table) {}
    
    virtual ~OperationEvalHessianBsplineClenshawCurtis() {}
    
    virtual double evalHessian(DataVector &alpha, const std::vector<double> &point,
                               DataVector &gradient, DataMatrix &hessian);
    
protected:
    GridStorage *storage;
    SBsplineClenshawCurtisBase base;
};

}
}

#endif /* OPERATIONEVALHESSIANBSPLINECLENSHAWCURTIS_HPP */
