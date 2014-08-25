/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_BASIS_NOBOUNDARY_OPERATION_OPERATIONEVALGRADIENTBSPLINE_HPP
#define SGPP_OPT_BASIS_NOBOUNDARY_OPERATION_OPERATIONEVALGRADIENTBSPLINE_HPP

#include "opt/operation/OperationEvalGradient.hpp"
#include "base/grid/GridStorage.hpp"
#include "opt/basis/bspline/noboundary/BsplineBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace opt
{

/**
 * Operation for evaluating B-spline linear combinations on Noboundary grids and their gradients.
 */
class OperationEvalGradientBspline : public OperationEvalGradient
{
public:
    /**
     * Constructor.
     * 
     * @param storage   storage of the sparse grid
     * @param degree    B-spline degree
     */
    OperationEvalGradientBspline(base::GridStorage *storage, size_t degree) :
        storage(storage), base(degree)
    {
    }
    
    /**
     * Virtual destructor.
     */
    virtual ~OperationEvalGradientBspline()
    {
    }
    
    /**
     * @param       alpha       coefficient vector
     * @param       point       evaluation point
     * @param[out]  gradient    gradient of linear combination
     * @return                  value of linear combination
     */
    virtual double evalGradient(base::DataVector &alpha, const std::vector<double> &point,
                                base::DataVector &gradient);
    
protected:
    /// storage of the sparse grid
    base::GridStorage *storage;
    /// 1D B-spline basis
    SBsplineBase base;
};

}
}

#endif
