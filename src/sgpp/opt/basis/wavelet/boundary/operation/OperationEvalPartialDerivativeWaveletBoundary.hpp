/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_BASIS_WAVELETBOUNDARY_NOBOUNDARY_OPERATION_OPERATIONEVALPARTIALDERIVATIVEWAVELETBOUNDARY_HPP
#define SGPP_OPT_BASIS_WAVELETBOUNDARY_NOBOUNDARY_OPERATION_OPERATIONEVALPARTIALDERIVATIVEWAVELETBOUNDARY_HPP

#include "opt/operation/OperationEvalPartialDerivative.hpp"
#include "base/grid/GridStorage.hpp"
#include "opt/basis/wavelet/boundary/WaveletBoundaryBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace opt
{

/**
 * Operation for evaluating partial derivatives of wavelet
 * linear combinations on Boundary grids.
 */
class OperationEvalPartialDerivativeWaveletBoundary :
        public OperationEvalPartialDerivative
{
public:
    /**
     * Constructor.
     * 
     * @param storage   storage of the sparse grid
     */
    OperationEvalPartialDerivativeWaveletBoundary(base::GridStorage *storage) :
        storage(storage)
    {
    }
    
    /**
     * Virtual destructor.
     */
    virtual ~OperationEvalPartialDerivativeWaveletBoundary()
    {
    }
    
    /**
     * @param alpha     coefficient vector
     * @param point     evaluation point
     * @param deriv_dim dimension in which the partial derivative should be taken
     * @return          value of the partial derivative of the linear combination
     */
    virtual double evalPartialDerivative(base::DataVector &alpha,
                                         const std::vector<double> &point,
                                         size_t deriv_dim);
    
protected:
    /// storage of the sparse grid
    base::GridStorage *storage;
    /// 1D wavelet basis
    SWaveletBoundaryBase base;
};

}
}

#endif
