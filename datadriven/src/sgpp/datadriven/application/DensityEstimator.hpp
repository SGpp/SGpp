// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DENSITYESTIMATOR_HPP_
#define DENSITYESTIMATOR_HPP_

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <iostream>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

class DensityEstimator {
public:
    DensityEstimator();
    DensityEstimator(base::DataMatrix& samples);
    virtual ~DensityEstimator();

    virtual void initialize(base::DataMatrix& samples) = 0;

    virtual float_t pdf(base::DataVector& x) = 0;
    virtual void pdf(base::DataMatrix& points, base::DataVector& res) = 0;

    virtual float_t mean() = 0;
    virtual float_t variance() = 0;
    virtual float_t std_deviation() = 0;
    virtual void cov(base::DataMatrix& cov) = 0;
    virtual void corrcoef(base::DataMatrix& corr);

    virtual base::DataVector* getSamples(size_t dim) = 0;
    virtual base::DataMatrix* getSamples();

    virtual size_t getDim() = 0;
    virtual size_t getNsamples() = 0;

protected:
    base::DataMatrix samples;
};

} /* namespace datadriven */
} /* namespace sg */

#endif /* DENSITYESTIMATOR_HPP_ */
