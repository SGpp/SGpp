/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_BASIS_BSPLINE_MODIFIED_OPERATION_OPERATIONEVALMODBSPLINE_HPP
#define SGPP_OPT_BASIS_BSPLINE_MODIFIED_OPERATION_OPERATIONEVALMODBSPLINE_HPP

#include "base/operation/OperationEval.hpp"
#include "base/grid/GridStorage.hpp"
#include "opt/basis/bspline/modified/ModBsplineBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace opt
{

/**
 * Operation for evaluating modified B-spline linear combinations on Noboundary grids.
 */
class OperationEvalModBspline : public base::OperationEval
{
public:
    /**
     * Constructor.
     * 
     * @param storage   storage of the sparse grid
     * @param degree    B-spline degree
     */
    OperationEvalModBspline(base::GridStorage *storage, size_t degree) :
        storage(storage), base(degree)
    {
    }
    
    /**
     * Virtual destructor.
     */
    virtual ~OperationEvalModBspline()
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
    SModBsplineBase base;
};

}
}

#endif
