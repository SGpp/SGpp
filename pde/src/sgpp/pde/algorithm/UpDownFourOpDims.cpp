// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/algorithm/UpDownFourOpDims.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

UpDownFourOpDims::UpDownFourOpDims(sgpp::base::GridStorage* storage, double***** coef)
    : storage(storage),
      coefs((*coef)),
      algoDims(storage->getAlgorithmicDimensions()),
      numAlgoDims_(storage->getAlgorithmicDimensions().size()) {
  generateMap();
}

UpDownFourOpDims::UpDownFourOpDims(sgpp::base::GridStorage* storage)
    : storage(storage),
      coefs(nullptr),
      algoDims(storage->getAlgorithmicDimensions()),
      numAlgoDims_(storage->getAlgorithmicDimensions().size()) {
  generateMap();
}

UpDownFourOpDims::~UpDownFourOpDims() {}

void UpDownFourOpDims::generateMap() {
  // Build the function mapping

  // unidirectional...
  fnMap.insert(std::make_pair(0, &sgpp::pde::UpDownFourOpDims::specialOpUnidirectional));

  // singles...
  fnMap.insert(std::make_pair(1, &sgpp::pde::UpDownFourOpDims::specialOpFour));
  fnMap.insert(std::make_pair(2, &sgpp::pde::UpDownFourOpDims::specialOpThree));
  fnMap.insert(std::make_pair(4, &sgpp::pde::UpDownFourOpDims::specialOpTwo));
  fnMap.insert(std::make_pair(8, &sgpp::pde::UpDownFourOpDims::specialOpOne));

  // doubles
  fnMap.insert(std::make_pair(3, &sgpp::pde::UpDownFourOpDims::specialOpThreeAndOpFour));
  fnMap.insert(std::make_pair(5, &sgpp::pde::UpDownFourOpDims::specialOpTwoAndOpFour));
  fnMap.insert(std::make_pair(6, &sgpp::pde::UpDownFourOpDims::specialOpTwoAndOpThree));
  fnMap.insert(std::make_pair(9, &sgpp::pde::UpDownFourOpDims::specialOpOneAndOpFour));
  fnMap.insert(std::make_pair(10, &sgpp::pde::UpDownFourOpDims::specialOpOneAndOpThree));
  fnMap.insert(std::make_pair(12, &sgpp::pde::UpDownFourOpDims::specialOpOneAndOpTwo));

  // triples
  fnMap.insert(std::make_pair(7, &sgpp::pde::UpDownFourOpDims::specialOpTwoAndOpThreeAndOpFour));
  fnMap.insert(std::make_pair(11, &sgpp::pde::UpDownFourOpDims::specialOpOneAndOpThreeAndOpFour));
  fnMap.insert(std::make_pair(13, &sgpp::pde::UpDownFourOpDims::specialOpOneAndOpTwoAndOpFour));
  fnMap.insert(std::make_pair(14, &sgpp::pde::UpDownFourOpDims::specialOpOneAndOpTwoAndOpThree));

  // quadruple
  fnMap.insert(
      std::make_pair(15, &sgpp::pde::UpDownFourOpDims::specialOpOneAndOpTwoAndOpThreeAndOpFour));
}

void UpDownFourOpDims::updown(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                              size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
                              size_t op_dim_four) {
  size_t num = 0;

  if (dim == op_dim_one) num += 8;

  if (dim == op_dim_two) num += 4;

  if (dim == op_dim_three) num += 2;

  if (dim == op_dim_four) num += 1;

  // Call the relevant function...
  MFP fp = fnMap[num];
  (this->*fp)(alpha, result, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  result.setAll(0.0);

#pragma omp parallel
  {
#pragma omp single nowait
    {
      for (size_t i = 0; i < this->numAlgoDims_; i++) {
        for (size_t j = 0; j < this->numAlgoDims_; j++) {
          for (size_t k = 0; k < this->numAlgoDims_; k++) {
            for (size_t l = 0; l < this->numAlgoDims_; l++) {
#pragma omp task firstprivate(i, j) shared(alpha, result)
              {
                sgpp::base::DataVector beta(result.getSize());

                if (this->coefs != nullptr) {
                  if (this->coefs[i][j][k][l] != 0.0) {
                    this->updown(alpha, beta, this->numAlgoDims_ - 1, i, j, k, l);

#pragma omp critical
                    { result.axpy(this->coefs[i][j][k][l], beta); }
                  }
                } else {
                  this->updown(alpha, beta, this->numAlgoDims_ - 1, i, j, k, l);

#pragma omp critical
                  { result.add(beta); }
                }
              }
            }
          }
        }
      }

#pragma omp taskwait
    }
  }
}

void UpDownFourOpDims::specialOpX(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
    void (sgpp::pde::UpDownFourOpDims::*pt2UpFunc)(sgpp::base::DataVector&, sgpp::base::DataVector&,
                                                   size_t),
    void (sgpp::pde::UpDownFourOpDims::*pt2DownFunc)(sgpp::base::DataVector&,
                                                     sgpp::base::DataVector&, size_t),
    size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
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
      (this->*pt2UpFunc)(alpha, temp, this->algoDims[dim]);
      updown(temp, result, dim - 1, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

// Same from the other direction:
#pragma omp task if (curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp_two, \
                                                                        result_temp)
    {  // NOLINT(whitespace/braces)
      updown(alpha, temp_two, dim - 1, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
      (this->*pt2DownFunc)(temp_two, result_temp, this->algoDims[dim]);
    }

#pragma omp taskwait

    result.add(result_temp);
  } else {
    // Terminates dimension recursion
    sgpp::base::DataVector temp(alpha.getSize());

#pragma omp task if (curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, result)
    (this->*pt2UpFunc)(alpha, result, this->algoDims[dim]);

#pragma omp task if (curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp)
    (this->*pt2DownFunc)(alpha, temp, this->algoDims[dim]);

#pragma omp taskwait

    result.add(temp);
  }
}

void UpDownFourOpDims::specialOpUnidirectional(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataVector& result, size_t dim,
                                               size_t op_dim_one, size_t op_dim_two,
                                               size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &sgpp::pde::UpDownFourOpDims::up, &sgpp::pde::UpDownFourOpDims::down,
             dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpOne(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                    size_t dim, size_t op_dim_one, size_t op_dim_two,
                                    size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &sgpp::pde::UpDownFourOpDims::upOpDimOne,
             &sgpp::pde::UpDownFourOpDims::downOpDimOne, dim, op_dim_one, op_dim_two, op_dim_three,
             op_dim_four);
}

void UpDownFourOpDims::specialOpTwo(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                    size_t dim, size_t op_dim_one, size_t op_dim_two,
                                    size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &sgpp::pde::UpDownFourOpDims::upOpDimTwo,
             &sgpp::pde::UpDownFourOpDims::downOpDimTwo, dim, op_dim_one, op_dim_two, op_dim_three,
             op_dim_four);
}

void UpDownFourOpDims::specialOpThree(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                      size_t dim, size_t op_dim_one, size_t op_dim_two,
                                      size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &sgpp::pde::UpDownFourOpDims::upOpDimThree,
             &sgpp::pde::UpDownFourOpDims::downOpDimThree, dim, op_dim_one, op_dim_two,
             op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpFour(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                     size_t dim, size_t op_dim_one, size_t op_dim_two,
                                     size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &sgpp::pde::UpDownFourOpDims::upOpDimFour,
             &sgpp::pde::UpDownFourOpDims::downOpDimFour, dim, op_dim_one, op_dim_two, op_dim_three,
             op_dim_four);
}

void UpDownFourOpDims::specialOpOneAndOpTwo(sgpp::base::DataVector& alpha,
                                            sgpp::base::DataVector& result, size_t dim,
                                            size_t op_dim_one, size_t op_dim_two,
                                            size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &sgpp::pde::UpDownFourOpDims::upOpDimOneAndOpDimTwo,
             &sgpp::pde::UpDownFourOpDims::downOpDimOneAndOpDimTwo, dim, op_dim_one, op_dim_two,
             op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpOneAndOpThree(sgpp::base::DataVector& alpha,
                                              sgpp::base::DataVector& result, size_t dim,
                                              size_t op_dim_one, size_t op_dim_two,
                                              size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &sgpp::pde::UpDownFourOpDims::upOpDimOneAndOpDimThree,
             &sgpp::pde::UpDownFourOpDims::downOpDimOneAndOpDimThree, dim, op_dim_one, op_dim_two,
             op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpOneAndOpFour(sgpp::base::DataVector& alpha,
                                             sgpp::base::DataVector& result, size_t dim,
                                             size_t op_dim_one, size_t op_dim_two,
                                             size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &sgpp::pde::UpDownFourOpDims::upOpDimOneAndOpDimFour,
             &sgpp::pde::UpDownFourOpDims::downOpDimOneAndOpDimFour, dim, op_dim_one, op_dim_two,
             op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpTwoAndOpThree(sgpp::base::DataVector& alpha,
                                              sgpp::base::DataVector& result, size_t dim,
                                              size_t op_dim_one, size_t op_dim_two,
                                              size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &sgpp::pde::UpDownFourOpDims::upOpDimTwoAndOpDimThree,
             &sgpp::pde::UpDownFourOpDims::downOpDimTwoAndOpDimThree, dim, op_dim_one, op_dim_two,
             op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpTwoAndOpFour(sgpp::base::DataVector& alpha,
                                             sgpp::base::DataVector& result, size_t dim,
                                             size_t op_dim_one, size_t op_dim_two,
                                             size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &sgpp::pde::UpDownFourOpDims::upOpDimTwoAndOpDimFour,
             &sgpp::pde::UpDownFourOpDims::downOpDimTwoAndOpDimFour, dim, op_dim_one, op_dim_two,
             op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpThreeAndOpFour(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataVector& result, size_t dim,
                                               size_t op_dim_one, size_t op_dim_two,
                                               size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &sgpp::pde::UpDownFourOpDims::upOpDimThreeAndOpDimFour,
             &sgpp::pde::UpDownFourOpDims::downOpDimThreeAndOpDimFour, dim, op_dim_one, op_dim_two,
             op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpOneAndOpTwoAndOpThree(sgpp::base::DataVector& alpha,
                                                      sgpp::base::DataVector& result, size_t dim,
                                                      size_t op_dim_one, size_t op_dim_two,
                                                      size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &sgpp::pde::UpDownFourOpDims::upOpDimOneAndOpDimTwoAndOpDimThree,
             &sgpp::pde::UpDownFourOpDims::downOpDimOneAndOpDimTwoAndOpDimThree, dim, op_dim_one,
             op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpOneAndOpTwoAndOpFour(sgpp::base::DataVector& alpha,
                                                     sgpp::base::DataVector& result, size_t dim,
                                                     size_t op_dim_one, size_t op_dim_two,
                                                     size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &sgpp::pde::UpDownFourOpDims::upOpDimOneAndOpDimTwoAndOpDimFour,
             &sgpp::pde::UpDownFourOpDims::downOpDimOneAndOpDimTwoAndOpDimFour, dim, op_dim_one,
             op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpOneAndOpThreeAndOpFour(sgpp::base::DataVector& alpha,
                                                       sgpp::base::DataVector& result, size_t dim,
                                                       size_t op_dim_one, size_t op_dim_two,
                                                       size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &sgpp::pde::UpDownFourOpDims::upOpDimOneAndOpDimThreeAndOpDimFour,
             &sgpp::pde::UpDownFourOpDims::downOpDimOneAndOpDimThreeAndOpDimFour, dim, op_dim_one,
             op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpTwoAndOpThreeAndOpFour(sgpp::base::DataVector& alpha,
                                                       sgpp::base::DataVector& result, size_t dim,
                                                       size_t op_dim_one, size_t op_dim_two,
                                                       size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &sgpp::pde::UpDownFourOpDims::upOpDimTwoAndOpDimThreeAndOpDimFour,
             &sgpp::pde::UpDownFourOpDims::downOpDimTwoAndOpDimThreeAndOpDimFour, dim, op_dim_one,
             op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpOneAndOpTwoAndOpThreeAndOpFour(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim, size_t op_dim_one,
    size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result,
             &sgpp::pde::UpDownFourOpDims::upOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour,
             &sgpp::pde::UpDownFourOpDims::downOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour, dim,
             op_dim_one, op_dim_two, op_dim_three, op_dim_four);
}
}  // namespace pde
}  // namespace sgpp
