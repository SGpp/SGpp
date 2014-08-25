/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_BASIS_LINEAR_MODIFIED_OPERATION_OPERATIONEVALMODLINEAR_HPP
#define SGPP_OPT_BASIS_LINEAR_MODIFIED_OPERATION_OPERATIONEVALMODLINEAR_HPP

#include "base/operation/OperationEval.hpp"
#include "base/grid/GridStorage.hpp"
#include "opt/basis/linear/modified/ModLinearBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace opt
{

class OperationEvalModLinear : public base::OperationEval
{
public:
    /**
     * Constructor.
     * 
     * @param storage   storage of the sparse grid
     */
    OperationEvalModLinear(base::GridStorage *storage) : storage(storage)
    {
    }
    
    /**
     * Virtual destructor.
     */
    virtual ~OperationEvalModLinear()
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
    /// 1D linear basis
    SModLinearBase base;
};

}
}

#endif
