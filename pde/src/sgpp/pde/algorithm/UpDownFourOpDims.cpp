// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/algorithm/UpDownFourOpDims.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace pde {

UpDownFourOpDims::UpDownFourOpDims(SGPP::base::GridStorage* storage,
                                   float_t**** * coef) : storage(storage), coefs((*coef)),
  algoDims(storage->getAlgorithmicDimensions()),
  numAlgoDims_(storage->getAlgorithmicDimensions().size()) {
  generateMap();
}

UpDownFourOpDims::UpDownFourOpDims(SGPP::base::GridStorage* storage) : storage(
    storage), coefs(NULL), algoDims(storage->getAlgorithmicDimensions()),
  numAlgoDims_(storage->getAlgorithmicDimensions().size()) {
  generateMap();
}

UpDownFourOpDims::~UpDownFourOpDims() {
}

void UpDownFourOpDims::generateMap() {
  // Build the function mapping

  // unidirectional...
  fnMap.insert( std::make_pair( 0,
                                &SGPP::pde::UpDownFourOpDims::specialOpUnidirectional ));

  // singles...
  fnMap.insert( std::make_pair( 1, &SGPP::pde::UpDownFourOpDims::specialOpFour ));
  fnMap.insert( std::make_pair( 2,
                                &SGPP::pde::UpDownFourOpDims::specialOpThree ));
  fnMap.insert( std::make_pair( 4, &SGPP::pde::UpDownFourOpDims::specialOpTwo ));
  fnMap.insert( std::make_pair( 8, &SGPP::pde::UpDownFourOpDims::specialOpOne ));

  // float_ts
  fnMap.insert( std::make_pair( 3,
                                &SGPP::pde::UpDownFourOpDims::specialOpThreeAndOpFour ));
  fnMap.insert( std::make_pair( 5,
                                &SGPP::pde::UpDownFourOpDims::specialOpTwoAndOpFour ));
  fnMap.insert( std::make_pair( 6,
                                &SGPP::pde::UpDownFourOpDims::specialOpTwoAndOpThree ));
  fnMap.insert( std::make_pair( 9,
                                &SGPP::pde::UpDownFourOpDims::specialOpOneAndOpFour ));
  fnMap.insert( std::make_pair( 10,
                                &SGPP::pde::UpDownFourOpDims::specialOpOneAndOpThree ));
  fnMap.insert( std::make_pair( 12,
                                &SGPP::pde::UpDownFourOpDims::specialOpOneAndOpTwo ));

  // triples
  fnMap.insert( std::make_pair( 7,
                                &SGPP::pde::UpDownFourOpDims::specialOpTwoAndOpThreeAndOpFour ));
  fnMap.insert( std::make_pair( 11,
                                &SGPP::pde::UpDownFourOpDims::specialOpOneAndOpThreeAndOpFour ));
  fnMap.insert( std::make_pair( 13,
                                &SGPP::pde::UpDownFourOpDims::specialOpOneAndOpTwoAndOpFour ));
  fnMap.insert( std::make_pair( 14,
                                &SGPP::pde::UpDownFourOpDims::specialOpOneAndOpTwoAndOpThree ));

  // quadruple
  fnMap.insert( std::make_pair( 15,
                                &SGPP::pde::UpDownFourOpDims::specialOpOneAndOpTwoAndOpThreeAndOpFour ));
}

void UpDownFourOpDims::updown(SGPP::base::DataVector& alpha,
                              SGPP::base::DataVector& result, size_t dim, size_t op_dim_one,
                              size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
  size_t num = 0;

  if (dim == op_dim_one)
    num += 8;

  if (dim == op_dim_two)
    num += 4;

  if (dim == op_dim_three)
    num += 2;

  if (dim == op_dim_four)
    num += 1;

  // Call the relevant function...
  MFP fp = fnMap[num];
  (this->*fp)(alpha, result, dim, op_dim_one, op_dim_two, op_dim_three,
              op_dim_four);
}

void UpDownFourOpDims::mult(SGPP::base::DataVector& alpha,
                            SGPP::base::DataVector& result) {
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
                SGPP::base::DataVector beta(result.getSize());

                if (this->coefs != NULL) {
                  if (this->coefs[i][j][k][l] != 0.0) {
                    this->updown(alpha, beta, this->numAlgoDims_ - 1, i, j, k, l);

                    #pragma omp critical
                    {
                      result.axpy(this->coefs[i][j][k][l], beta);
                    }
                  }
                } else
                {
                  this->updown(alpha, beta, this->numAlgoDims_ - 1, i, j, k, l);

                  #pragma omp critical
                  {
                    result.add(beta);
                  }
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

void UpDownFourOpDims::specialOpX(SGPP::base::DataVector& alpha,
                                  SGPP::base::DataVector& result,
                                  void (SGPP::pde::UpDownFourOpDims::*pt2UpFunc)(SGPP::base::DataVector&,
                                      SGPP::base::DataVector&, size_t),
                                  void (SGPP::pde::UpDownFourOpDims::*pt2DownFunc)(SGPP::base::DataVector&,
                                      SGPP::base::DataVector&, size_t), size_t dim, size_t op_dim_one,
                                  size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
  size_t curNumAlgoDims = this->numAlgoDims_;
  size_t curMaxParallelDims = this->maxParallelDims_;

  //Unidirectional scheme
  if (dim > 0) {
    // Reordering ups and downs
    SGPP::base::DataVector temp(alpha.getSize());
    SGPP::base::DataVector result_temp(alpha.getSize());
    SGPP::base::DataVector temp_two(alpha.getSize());

    #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp, result)
    {
      (this->*pt2UpFunc)(alpha, temp, this->algoDims[dim]);
      updown(temp, result, dim - 1, op_dim_one, op_dim_two, op_dim_three,
             op_dim_four);
    }

    // Same from the other direction:
    #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp_two, result_temp)
    {
      updown(alpha, temp_two, dim - 1, op_dim_one, op_dim_two, op_dim_three,
             op_dim_four);
      (this->*pt2DownFunc)(temp_two, result_temp, this->algoDims[dim]);
    }

    #pragma omp taskwait

    result.add(result_temp);
  } else {
    // Terminates dimension recursion
    SGPP::base::DataVector temp(alpha.getSize());

    #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, result)
    (this->*pt2UpFunc)(alpha, result, this->algoDims[dim]);

    #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp)
    (this->*pt2DownFunc)(alpha, temp, this->algoDims[dim]);

    #pragma omp taskwait

    result.add(temp);
  }
}

void UpDownFourOpDims::specialOpUnidirectional(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result, size_t dim, size_t op_dim_one,
    size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &SGPP::pde::UpDownFourOpDims::up,
             &SGPP::pde::UpDownFourOpDims::down, dim, op_dim_one, op_dim_two, op_dim_three,
             op_dim_four);
}

void UpDownFourOpDims::specialOpOne(SGPP::base::DataVector& alpha,
                                    SGPP::base::DataVector& result, size_t dim, size_t op_dim_one,
                                    size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &SGPP::pde::UpDownFourOpDims::upOpDimOne,
             &SGPP::pde::UpDownFourOpDims::downOpDimOne, dim, op_dim_one, op_dim_two,
             op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpTwo(SGPP::base::DataVector& alpha,
                                    SGPP::base::DataVector& result, size_t dim, size_t op_dim_one,
                                    size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &SGPP::pde::UpDownFourOpDims::upOpDimTwo,
             &SGPP::pde::UpDownFourOpDims::downOpDimTwo, dim, op_dim_one, op_dim_two,
             op_dim_three, op_dim_four);

}

void UpDownFourOpDims::specialOpThree(SGPP::base::DataVector& alpha,
                                      SGPP::base::DataVector& result, size_t dim, size_t op_dim_one,
                                      size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &SGPP::pde::UpDownFourOpDims::upOpDimThree,
             &SGPP::pde::UpDownFourOpDims::downOpDimThree, dim, op_dim_one, op_dim_two,
             op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpFour(SGPP::base::DataVector& alpha,
                                     SGPP::base::DataVector& result, size_t dim, size_t op_dim_one,
                                     size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &SGPP::pde::UpDownFourOpDims::upOpDimFour,
             &SGPP::pde::UpDownFourOpDims::downOpDimFour, dim, op_dim_one, op_dim_two,
             op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpOneAndOpTwo(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result, size_t dim, size_t op_dim_one,
    size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &SGPP::pde::UpDownFourOpDims::upOpDimOneAndOpDimTwo,
             &SGPP::pde::UpDownFourOpDims::downOpDimOneAndOpDimTwo, dim, op_dim_one,
             op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpOneAndOpThree(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result, size_t dim, size_t op_dim_one,
    size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &SGPP::pde::UpDownFourOpDims::upOpDimOneAndOpDimThree,
             &SGPP::pde::UpDownFourOpDims::downOpDimOneAndOpDimThree, dim, op_dim_one,
             op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpOneAndOpFour(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result, size_t dim, size_t op_dim_one,
    size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &SGPP::pde::UpDownFourOpDims::upOpDimOneAndOpDimFour,
             &SGPP::pde::UpDownFourOpDims::downOpDimOneAndOpDimFour, dim, op_dim_one,
             op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpTwoAndOpThree(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result, size_t dim, size_t op_dim_one,
    size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &SGPP::pde::UpDownFourOpDims::upOpDimTwoAndOpDimThree,
             &SGPP::pde::UpDownFourOpDims::downOpDimTwoAndOpDimThree, dim, op_dim_one,
             op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpTwoAndOpFour(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result, size_t dim, size_t op_dim_one,
    size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result, &SGPP::pde::UpDownFourOpDims::upOpDimTwoAndOpDimFour,
             &SGPP::pde::UpDownFourOpDims::downOpDimTwoAndOpDimFour, dim, op_dim_one,
             op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpThreeAndOpFour(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result, size_t dim, size_t op_dim_one,
    size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
  specialOpX(alpha, result,
             &SGPP::pde::UpDownFourOpDims::upOpDimThreeAndOpDimFour,
             &SGPP::pde::UpDownFourOpDims::downOpDimThreeAndOpDimFour, dim, op_dim_one,
             op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpOneAndOpTwoAndOpThree(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim,
  size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
  size_t op_dim_four) {
  specialOpX(alpha, result,
             &SGPP::pde::UpDownFourOpDims::upOpDimOneAndOpDimTwoAndOpDimThree,
             &SGPP::pde::UpDownFourOpDims::downOpDimOneAndOpDimTwoAndOpDimThree, dim,
             op_dim_one, op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpOneAndOpTwoAndOpFour(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim,
  size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
  size_t op_dim_four) {
  specialOpX(alpha, result,
             &SGPP::pde::UpDownFourOpDims::upOpDimOneAndOpDimTwoAndOpDimFour,
             &SGPP::pde::UpDownFourOpDims::downOpDimOneAndOpDimTwoAndOpDimFour, dim,
             op_dim_one, op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpOneAndOpThreeAndOpFour(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim,
  size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
  size_t op_dim_four) {
  specialOpX(alpha, result,
             &SGPP::pde::UpDownFourOpDims::upOpDimOneAndOpDimThreeAndOpDimFour,
             &SGPP::pde::UpDownFourOpDims::downOpDimOneAndOpDimThreeAndOpDimFour, dim,
             op_dim_one, op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpTwoAndOpThreeAndOpFour(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim,
  size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
  size_t op_dim_four) {
  specialOpX(alpha, result,
             &SGPP::pde::UpDownFourOpDims::upOpDimTwoAndOpDimThreeAndOpDimFour,
             &SGPP::pde::UpDownFourOpDims::downOpDimTwoAndOpDimThreeAndOpDimFour, dim,
             op_dim_one, op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpOneAndOpTwoAndOpThreeAndOpFour(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim,
  size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
  size_t op_dim_four) {
  specialOpX(alpha, result,
             &SGPP::pde::UpDownFourOpDims::upOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour,
             &SGPP::pde::UpDownFourOpDims::downOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour,
             dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
}



}
}