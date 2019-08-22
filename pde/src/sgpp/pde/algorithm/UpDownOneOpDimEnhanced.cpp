// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/algorithm/UpDownOneOpDimEnhanced.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

UpDownOneOpDimEnhanced::UpDownOneOpDimEnhanced(sgpp::base::GridStorage* storage,
                                               sgpp::base::DataVector& coef)
    : storage(storage),
      coefs(&coef),
      algoDims(storage->getAlgorithmicDimensions()),
      numAlgoDims_(storage->getAlgorithmicDimensions().size()) {}

UpDownOneOpDimEnhanced::UpDownOneOpDimEnhanced(sgpp::base::GridStorage* storage)
    : storage(storage),
      coefs(nullptr),
      algoDims(storage->getAlgorithmicDimensions()),
      numAlgoDims_(storage->getAlgorithmicDimensions().size()) {}

UpDownOneOpDimEnhanced::~UpDownOneOpDimEnhanced() {}

void UpDownOneOpDimEnhanced::mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  sgpp::base::DataMatrix beta(result.getSize(), this->numAlgoDims_);
  sgpp::base::DataMatrix maAlpha(alpha.getSize(), this->numAlgoDims_);

  result.setAll(0.0);
  maAlpha.expand(alpha);

#pragma omp parallel
  {
#pragma omp single nowait
    { this->updown(maAlpha, beta, this->numAlgoDims_ - 1); }
  }

  if (coefs == nullptr) {
    beta.addReduce(result);
  } else {
    beta.addReduce(result, *coefs, 0);
  }
}

void UpDownOneOpDimEnhanced::multParallelBuildingBlock(sgpp::base::DataVector& alpha,
                                                       sgpp::base::DataVector& result) {
  sgpp::base::DataMatrix beta(result.getSize(), this->numAlgoDims_);
  sgpp::base::DataMatrix maAlpha(alpha.getSize(), this->numAlgoDims_);

  result.setAll(0.0);
  maAlpha.expand(alpha);

  this->updown(maAlpha, beta, this->numAlgoDims_ - 1);

  if (coefs == nullptr) {
    beta.addReduce(result);
  } else {
    beta.addReduce(result, *coefs, 0);
  }
}

void UpDownOneOpDimEnhanced::updown(sgpp::base::DataMatrix& alpha, sgpp::base::DataMatrix& result,
                                    size_t dim) {
  size_t curNumAlgoDims = this->numAlgoDims_;
  size_t curMaxParallelDims = this->maxParallelDims_;

  // Unidirectional scheme
  if (dim > 0) {
    // Reordering ups and downs
    sgpp::base::DataMatrix temp(alpha.getNrows(), this->numAlgoDims_);
    sgpp::base::DataMatrix result_temp(alpha.getNrows(), this->numAlgoDims_);
    sgpp::base::DataMatrix temp_two(alpha.getNrows(), this->numAlgoDims_);

#pragma omp task if (curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp, result)
    {
      up(alpha, temp, dim);
      updown(temp, result, dim - 1);
    }

#pragma omp task if (curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp_two, \
                                                                        result_temp)
    {  // NOLINT(whitespace/braces)
      updown(alpha, temp_two, dim - 1);
      down(temp_two, result_temp, dim);
    }

#pragma omp taskwait

    result.add(result_temp);
  } else {
    // Terminates dimension recursion
    sgpp::base::DataMatrix temp(alpha.getNrows(), this->numAlgoDims_);

#pragma omp task if (curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, result)
    up(alpha, result, dim);

#pragma omp task if (curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp)
    down(alpha, temp, dim);

#pragma omp taskwait

    result.add(temp);
  }
}
}  // namespace pde
}  // namespace sgpp
