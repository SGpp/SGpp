/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Sam Maurus (MA thesis)

#include <sgpp/pde/algorithm/UpDownFourOpDims.hpp>

namespace sg {
  namespace pde {

    UpDownFourOpDims::UpDownFourOpDims(sg::base::GridStorage* storage, double**** * coef) : storage(storage), coefs((*coef)), algoDims(storage->getAlgorithmicDimensions()), numAlgoDims_(storage->getAlgorithmicDimensions().size()) {
      generateMap();
    }

    UpDownFourOpDims::UpDownFourOpDims(sg::base::GridStorage* storage) : storage(storage), coefs(NULL), algoDims(storage->getAlgorithmicDimensions()), numAlgoDims_(storage->getAlgorithmicDimensions().size()) {
      generateMap();
    }

    UpDownFourOpDims::~UpDownFourOpDims() {
    }

    void UpDownFourOpDims::generateMap() {
      // Build the function mapping

      // unidirectional...
      fnMap.insert( std::make_pair( 0, &sg::pde::UpDownFourOpDims::specialOpUnidirectional ));

      // singles...
      fnMap.insert( std::make_pair( 1, &sg::pde::UpDownFourOpDims::specialOpFour ));
      fnMap.insert( std::make_pair( 2, &sg::pde::UpDownFourOpDims::specialOpThree ));
      fnMap.insert( std::make_pair( 4, &sg::pde::UpDownFourOpDims::specialOpTwo ));
      fnMap.insert( std::make_pair( 8, &sg::pde::UpDownFourOpDims::specialOpOne ));

      // doubles
      fnMap.insert( std::make_pair( 3, &sg::pde::UpDownFourOpDims::specialOpThreeAndOpFour ));
      fnMap.insert( std::make_pair( 5, &sg::pde::UpDownFourOpDims::specialOpTwoAndOpFour ));
      fnMap.insert( std::make_pair( 6, &sg::pde::UpDownFourOpDims::specialOpTwoAndOpThree ));
      fnMap.insert( std::make_pair( 9, &sg::pde::UpDownFourOpDims::specialOpOneAndOpFour ));
      fnMap.insert( std::make_pair( 10, &sg::pde::UpDownFourOpDims::specialOpOneAndOpThree ));
      fnMap.insert( std::make_pair( 12, &sg::pde::UpDownFourOpDims::specialOpOneAndOpTwo ));

      // triples
      fnMap.insert( std::make_pair( 7, &sg::pde::UpDownFourOpDims::specialOpTwoAndOpThreeAndOpFour ));
      fnMap.insert( std::make_pair( 11, &sg::pde::UpDownFourOpDims::specialOpOneAndOpThreeAndOpFour ));
      fnMap.insert( std::make_pair( 13, &sg::pde::UpDownFourOpDims::specialOpOneAndOpTwoAndOpFour ));
      fnMap.insert( std::make_pair( 14, &sg::pde::UpDownFourOpDims::specialOpOneAndOpTwoAndOpThree ));

      // quadruple
      fnMap.insert( std::make_pair( 15, &sg::pde::UpDownFourOpDims::specialOpOneAndOpTwoAndOpThreeAndOpFour ));
    }

    void UpDownFourOpDims::updown(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
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
      (this->*fp)(alpha, result, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

    void UpDownFourOpDims::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
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
                    sg::base::DataVector beta(result.getSize());

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

    void UpDownFourOpDims::specialOpX(sg::base::DataVector& alpha, sg::base::DataVector& result, void (sg::pde::UpDownFourOpDims::*pt2UpFunc)(sg::base::DataVector&, sg::base::DataVector&, size_t), void (sg::pde::UpDownFourOpDims::*pt2DownFunc)(sg::base::DataVector&, sg::base::DataVector&, size_t), size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
      size_t curNumAlgoDims = this->numAlgoDims_;
      size_t curMaxParallelDims = this->maxParallelDims_;

      //Unidirectional scheme
      if (dim > 0) {
        // Reordering ups and downs
        sg::base::DataVector temp(alpha.getSize());
        sg::base::DataVector result_temp(alpha.getSize());
        sg::base::DataVector temp_two(alpha.getSize());

        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp, result)
        {
          (this->*pt2UpFunc)(alpha, temp, this->algoDims[dim]);
          updown(temp, result, dim - 1, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
        }

        // Same from the other direction:
        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp_two, result_temp)
        {
          updown(alpha, temp_two, dim - 1, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
          (this->*pt2DownFunc)(temp_two, result_temp, this->algoDims[dim]);
        }

        #pragma omp taskwait

        result.add(result_temp);
      } else {
        // Terminates dimension recursion
        sg::base::DataVector temp(alpha.getSize());

        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, result)
        (this->*pt2UpFunc)(alpha, result, this->algoDims[dim]);

        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp)
        (this->*pt2DownFunc)(alpha, temp, this->algoDims[dim]);

        #pragma omp taskwait

        result.add(temp);
      }
    }

    void UpDownFourOpDims::specialOpUnidirectional(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::up, &sg::pde::UpDownFourOpDims::down, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

    void UpDownFourOpDims::specialOpOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimOne, &sg::pde::UpDownFourOpDims::downOpDimOne, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

    void UpDownFourOpDims::specialOpTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimTwo, &sg::pde::UpDownFourOpDims::downOpDimTwo, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);

    }

    void UpDownFourOpDims::specialOpThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimThree, &sg::pde::UpDownFourOpDims::downOpDimThree, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

    void UpDownFourOpDims::specialOpFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimFour, &sg::pde::UpDownFourOpDims::downOpDimFour, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

    void UpDownFourOpDims::specialOpOneAndOpTwo(sg::base::DataVector& alpha,
        sg::base::DataVector& result, size_t dim, size_t op_dim_one,
        size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimOneAndOpDimTwo, &sg::pde::UpDownFourOpDims::downOpDimOneAndOpDimTwo, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

    void UpDownFourOpDims::specialOpOneAndOpThree(sg::base::DataVector& alpha,
        sg::base::DataVector& result, size_t dim, size_t op_dim_one,
        size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimOneAndOpDimThree, &sg::pde::UpDownFourOpDims::downOpDimOneAndOpDimThree, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

    void UpDownFourOpDims::specialOpOneAndOpFour(sg::base::DataVector& alpha,
        sg::base::DataVector& result, size_t dim, size_t op_dim_one,
        size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimOneAndOpDimFour, &sg::pde::UpDownFourOpDims::downOpDimOneAndOpDimFour, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

    void UpDownFourOpDims::specialOpTwoAndOpThree(sg::base::DataVector& alpha,
        sg::base::DataVector& result, size_t dim, size_t op_dim_one,
        size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimTwoAndOpDimThree, &sg::pde::UpDownFourOpDims::downOpDimTwoAndOpDimThree, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

    void UpDownFourOpDims::specialOpTwoAndOpFour(sg::base::DataVector& alpha,
        sg::base::DataVector& result, size_t dim, size_t op_dim_one,
        size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimTwoAndOpDimFour, &sg::pde::UpDownFourOpDims::downOpDimTwoAndOpDimFour, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

    void UpDownFourOpDims::specialOpThreeAndOpFour(sg::base::DataVector& alpha,
        sg::base::DataVector& result, size_t dim, size_t op_dim_one,
        size_t op_dim_two, size_t op_dim_three, size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimThreeAndOpDimFour, &sg::pde::UpDownFourOpDims::downOpDimThreeAndOpDimFour, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

    void UpDownFourOpDims::specialOpOneAndOpTwoAndOpThree(
      sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim,
      size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
      size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimOneAndOpDimTwoAndOpDimThree, &sg::pde::UpDownFourOpDims::downOpDimOneAndOpDimTwoAndOpDimThree, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

    void UpDownFourOpDims::specialOpOneAndOpTwoAndOpFour(
      sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim,
      size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
      size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimOneAndOpDimTwoAndOpDimFour, &sg::pde::UpDownFourOpDims::downOpDimOneAndOpDimTwoAndOpDimFour, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

    void UpDownFourOpDims::specialOpOneAndOpThreeAndOpFour(
      sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim,
      size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
      size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimOneAndOpDimThreeAndOpDimFour, &sg::pde::UpDownFourOpDims::downOpDimOneAndOpDimThreeAndOpDimFour, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

    void UpDownFourOpDims::specialOpTwoAndOpThreeAndOpFour(
      sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim,
      size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
      size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimTwoAndOpDimThreeAndOpDimFour, &sg::pde::UpDownFourOpDims::downOpDimTwoAndOpDimThreeAndOpDimFour, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }

    void UpDownFourOpDims::specialOpOneAndOpTwoAndOpThreeAndOpFour(
      sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim,
      size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
      size_t op_dim_four) {
      specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour, &sg::pde::UpDownFourOpDims::downOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
    }



  }
}
