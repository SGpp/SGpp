/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de
// some defines for the following algorithm

#ifndef OPERATIONROSENBLATTTRANSFORMATIONKDE_HPP
#define OPERATIONROSENBLATTTRANSFORMATIONKDE_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/application/GaussianKDE.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

/**
 * Do transformation in all dimensions
 */
class OperationRosenblattTransformationKDE {
public:
    OperationRosenblattTransformationKDE(GaussianKDE& kde);
    virtual ~OperationRosenblattTransformationKDE();

    /**
     * Rosenblatt Transformation with mixed starting dimension
     *
     * @param pointsUniform data points to be transformed base::DataMatrix (rows: # of samples, columns: # of dims)
     * @param pointsCdf Output base::DataMatrix (rows: # of samples, columns: # of dims)
     */
    virtual void doTransformation(base::DataMatrix& pointsCdf,
            base::DataMatrix& pointsUniform);

//    /**
//     * Rosenblatt Transformation with fixed starting dimension
//     *
//     * @param pointsUniform data points to be transformed base::DataMatrix (rows: # of samples, columns: # of dims)
//     * @param pointsCdf Output base::DataMatrix (rows: # of samples, columns: # of dims)
//     */
//    virtual void doShuffeledTransformation(base::DataMatrix& pointsCdf,
//            base::DataMatrix& pointsUniform);

    /**
     *
     * @param x
     * @param samples1d
     * @param sigma
     * @param kern
     */
    float_t doTransformation1D(float_t x, base::DataVector& samples1d,
            float_t sigma, base::DataVector& kern);

private:
    datadriven::GaussianKDE* kde;
    base::DataVector bandwidths;

    size_t ndim;
    size_t nsamples;
};

}
}
#endif /* OPERATIONINVERSEROSENBLATTTRANSFORMATIONKDE_HPP */
