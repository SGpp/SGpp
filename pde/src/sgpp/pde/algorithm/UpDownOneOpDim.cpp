// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

UpDownOneOpDim::UpDownOneOpDim(sgpp::base::GridStorage* storage, sgpp::base::DataVector& coef)
    : storage(storage),
      coefs(&coef),
      algoDims(storage->getAlgorithmicDimensions()),
      numAlgoDims_(storage->getAlgorithmicDimensions().size()) {}

UpDownOneOpDim::UpDownOneOpDim(sgpp::base::GridStorage* storage)
    : storage(storage),
      coefs(nullptr),
      algoDims(storage->getAlgorithmicDimensions()),
      numAlgoDims_(storage->getAlgorithmicDimensions().size()) {}

UpDownOneOpDim::~UpDownOneOpDim() {}

void UpDownOneOpDim::mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  result.setAll(0.0);

#pragma omp parallel
  {
#pragma omp single nowait
    {
      for (size_t i = 0; i < this->numAlgoDims_; i++) {
#pragma omp task firstprivate(i) shared(alpha, result)
        {
          sgpp::base::DataVector beta(result.getSize());

          if (this->coefs != nullptr) {
            if (this->coefs->get(i) != 0.0) {
              this->updown(alpha, beta, this->numAlgoDims_ - 1, i);

#pragma omp critical
              { result.axpy(this->coefs->get(i), beta); }
            }
          } else {
            this->updown(alpha, beta, this->numAlgoDims_ - 1, i);

#pragma omp critical
            { result.add(beta); }
          }
        }
      }

#pragma omp taskwait
    }
  }
}

void UpDownOneOpDim::multParallelBuildingBlock(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataVector& result,
                                               size_t operationDim) {
  result.setAll(0.0);

  sgpp::base::DataVector beta(result.getSize());

  if (this->coefs != nullptr) {
    if (this->coefs->get(operationDim) != 0.0) {
      this->updown(alpha, beta, this->numAlgoDims_ - 1, operationDim);

      result.axpy(this->coefs->get(operationDim), beta);
    }
  } else {
    this->updown(alpha, beta, this->numAlgoDims_ - 1, operationDim);

    result.add(beta);
  }
}

void UpDownOneOpDim::updown(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                            size_t dim, size_t op_dim) {
  size_t curNumAlgoDims = this->numAlgoDims_;
  size_t curMaxParallelDims = this->maxParallelDims_;

  if (dim == op_dim) {
    specialOP(alpha, result, dim, op_dim);
  } else {
    // Unidirectional scheme
    if (dim > 0) {
      // Reordering ups and downs
      sgpp::base::DataVector temp(alpha.getSize());
      sgpp::base::DataVector result_temp(alpha.getSize());
      sgpp::base::DataVector temp_two(alpha.getSize());

#pragma omp task if (curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp, result)
      {
        up(alpha, temp, this->algoDims[dim]);
        updown(temp, result, dim - 1, op_dim);
      }

// Same from the other direction:
#pragma omp task if (curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp_two, \
                                                                        result_temp)
      {  // NOLINT(whitespace/braces)
        updown(alpha, temp_two, dim - 1, op_dim);
        down(temp_two, result_temp, this->algoDims[dim]);
      }

#pragma omp taskwait

      result.add(result_temp);
    } else {
      // Terminates dimension recursion
      sgpp::base::DataVector temp(alpha.getSize());

#pragma omp task if (curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, result)
      up(alpha, result, this->algoDims[dim]);

#pragma omp task if (curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp)
      down(alpha, temp, this->algoDims[dim]);

#pragma omp taskwait

      result.add(temp);
    }
  }
}

void UpDownOneOpDim::specialOP(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                               size_t dim, size_t op_dim) {
  size_t curNumAlgoDims = this->numAlgoDims_;
  size_t curMaxParallelDims = this->maxParallelDims_;

  // Unidirectional scheme
  if (dim > 0) {
    // Reordering ups and downs
    sgpp::base::DataVector temp(alpha.getSize());
    sgpp::base::DataVector result_temp(alpha.getSize());
    sgpp::base::DataVector temp_two(alpha.getSize());

#pragma omp task if (curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp, result)
    {
      upOpDim(alpha, temp, this->algoDims[dim]);
      updown(temp, result, dim - 1, op_dim);
    }

// Same from the other direction:
#pragma omp task if (curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp_two, \
                                                                        result_temp)
    {  // NOLINT(whitespace/braces)
      updown(alpha, temp_two, dim - 1, op_dim);
      downOpDim(temp_two, result_temp, this->algoDims[dim]);
    }

#pragma omp taskwait

    result.add(result_temp);
  } else {
    // Terminates dimension recursion
    sgpp::base::DataVector temp(alpha.getSize());

#pragma omp task if (curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, result)
    upOpDim(alpha, result, this->algoDims[dim]);

#pragma omp task if (curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp)
    downOpDim(alpha, temp, this->algoDims[dim]);

#pragma omp taskwait

    result.add(temp);
  }
}
}  // namespace pde
}  // namespace sgpp
