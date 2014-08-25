/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_BASIS_WAVELET_MODIFIED_OPERATION_OPERATIONEVALHESSIANMODWAVELET_HPP
#define SGPP_OPT_BASIS_WAVELET_MODIFIED_OPERATION_OPERATIONEVALHESSIANMODWAVELET_HPP

#include "opt/operation/OperationEvalHessian.hpp"
#include "base/grid/GridStorage.hpp"
#include "opt/basis/wavelet/modified/ModWaveletBasis.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

namespace sg
{
namespace opt
{

/**
 * Operation for evaluating modified wavelet linear combinations on Noboundary grids,
 * their gradients and their Hessians.
 */
class OperationEvalHessianModWavelet : public OperationEvalHessian
{
public:
    /**
     * Constructor.
     * 
     * @param storage   storage of the sparse grid
     */
    OperationEvalHessianModWavelet(base::GridStorage *storage) : storage(storage)
    {
    }
    
    /**
     * Virtual destructor.
     */
    virtual ~OperationEvalHessianModWavelet()
    {
    }
    
    /**
     * @param       alpha       coefficient vector
     * @param       point       evaluation point
     * @param[out]  gradient    gradient vector of linear combination
     * @param[out]  hessian     Hessian matrix of linear combination
     * @return                  value of linear combination
     */
    virtual double evalHessian(base::DataVector &alpha, const std::vector<double> &point,
                               base::DataVector &gradient, base::DataMatrix &hessian);
    
protected:
    /// storage of the sparse grid
    base::GridStorage *storage;
    /// 1D wavelet basis
    SModWaveletBase base;
};

}
}

#endif
