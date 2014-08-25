/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_BASIS_NOBOUNDARY_OPERATION_OPERATIONEVALBSPLINE_HPP
#define SGPP_OPT_BASIS_NOBOUNDARY_OPERATION_OPERATIONEVALBSPLINE_HPP

#include "base/operation/OperationEval.hpp"
#include "base/grid/GridStorage.hpp"
#include "opt/basis/bspline/noboundary/BsplineBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace opt
{

/**
 * Operation for evaluating B-spline linear combinations on Noboundary grids.
 */
class OperationEvalBspline : public base::OperationEval
{
public:
    /**
     * Constructor.
     * 
     * @param storage   storage of the sparse grid
     * @param degree    B-spline degree
     */
    OperationEvalBspline(base::GridStorage *storage, size_t degree) :
        storage(storage), base(degree)
    {
    }
    
    /**
     * Virtual destructor.
     */
    virtual ~OperationEvalBspline()
    {
    }
    
    /**
     * @param alpha     coefficient vector
     * @param point     evaluation point
     * @return          value of linear combination
     */
    virtual double eval(base::DataVector &alpha, std::vector<double> &point);
    
protected:
    /// storage of the sparse grid
    base::GridStorage *storage;
    /// 1D B-spline basis
    SBsplineBase base;
};

}
}

#endif
