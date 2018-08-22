// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {

class DensityEstimator {
 public:
  DensityEstimator();
  virtual ~DensityEstimator();

  virtual void initialize(base::DataMatrix& samples) = 0;

  virtual double pdf(base::DataVector& x) = 0;
  virtual void pdf(base::DataMatrix& points, base::DataVector& res) = 0;

  virtual double mean() = 0;
  virtual double variance() = 0;
  virtual double std_deviation();
  virtual void cov(base::DataMatrix& cov, base::DataMatrix* bounds = nullptr) = 0;
  virtual void corrcoef(base::DataMatrix& corr, base::DataMatrix* bounds = nullptr);

  virtual std::shared_ptr<base::DataVector> getSamples(size_t dim) = 0;
  virtual std::shared_ptr<base::DataMatrix> getSamples() = 0;

  virtual size_t getDim() = 0;
  virtual size_t getNsamples() = 0;

  double crossEntropy(sgpp::base::DataMatrix& samples);
};

}  // namespace datadriven
}  // namespace sgpp
