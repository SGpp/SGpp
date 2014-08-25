/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_BASIS_WAVELET_BOUNDARY_OPERATION_OPERATIONEVALGRADIENTWAVELETBOUNDARY_HPP
#define SGPP_OPT_BASIS_WAVELET_BOUNDARY_OPERATION_OPERATIONEVALGRADIENTWAVELETBOUNDARY_HPP

#include "opt/operation/OperationEvalGradient.hpp"
#include "base/grid/GridStorage.hpp"
#include "opt/basis/wavelet/boundary/WaveletBoundaryBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace opt
{

/**
 * Operation for evaluating wavelet linear combinations on Boundary grids and their gradients.
 */
class OperationEvalGradientWaveletBoundary : public OperationEvalGradient
{
public:
    /**
     * Constructor.
     * 
     * @param storage   storage of the sparse grid
     */
    OperationEvalGradientWaveletBoundary(base::GridStorage *storage) : storage(storage)
    {
    }
    
    /**
     * Virtual destructor.
     */
    virtual ~OperationEvalGradientWaveletBoundary()
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
    /// 1D wavelet basis
    SWaveletBoundaryBase base;
};

}
}

#endif
