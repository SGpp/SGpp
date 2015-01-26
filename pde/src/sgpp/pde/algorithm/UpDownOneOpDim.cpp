/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/algorithm/UpDownOneOpDim.hpp"

namespace sg {
  namespace pde {

    UpDownOneOpDim::UpDownOneOpDim(sg::base::GridStorage* storage, sg::base::DataVector& coef) : storage(storage), coefs(&coef), algoDims(storage->getAlgorithmicDimensions()), numAlgoDims_(storage->getAlgorithmicDimensions().size()) {
    }

    UpDownOneOpDim::UpDownOneOpDim(sg::base::GridStorage* storage): storage(storage), coefs(NULL), algoDims(storage->getAlgorithmicDimensions()), numAlgoDims_(storage->getAlgorithmicDimensions().size()) {
    }

    UpDownOneOpDim::~UpDownOneOpDim() {
    }

    void UpDownOneOpDim::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      result.setAll(0.0);

      #pragma omp parallel
      {
        #pragma omp single nowait
        {
          for (size_t i = 0; i < this->numAlgoDims_; i++) {
            #pragma omp task firstprivate(i) shared(alpha, result)
            {
              sg::base::DataVector beta(result.getSize());

              if (this->coefs != NULL) {
                if (this->coefs->get(i) != 0.0) {
                  this->updown(alpha, beta, this->numAlgoDims_ - 1, i);

                  #pragma omp critical
                  {
                    result.axpy(this->coefs->get(i), beta);
                  }
                }
              } else
              {
                this->updown(alpha, beta, this->numAlgoDims_ - 1, i);

                #pragma omp critical
                {
                  result.add(beta);
                }
              }
            }
          }

          #pragma omp taskwait
        }
      }

    }

    void UpDownOneOpDim::multParallelBuildingBlock(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t operationDim) {
      result.setAll(0.0);

      sg::base::DataVector beta(result.getSize());

      if (this->coefs != NULL) {
        if (this->coefs->get(operationDim) != 0.0) {
          this->updown(alpha, beta, this->numAlgoDims_ - 1, operationDim);

          result.axpy(this->coefs->get(operationDim), beta);
        }
      } else {
        this->updown(alpha, beta, this->numAlgoDims_ - 1, operationDim);

        result.add(beta);
      }
    }

    void UpDownOneOpDim::updown(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim) {
      size_t curNumAlgoDims = this->numAlgoDims_;
      size_t curMaxParallelDims = this->maxParallelDims_;

      if (dim == op_dim) {
        specialOP(alpha, result, dim, op_dim);
      } else {
        //Unidirectional scheme
        if (dim > 0) {
          // Reordering ups and downs
          sg::base::DataVector temp(alpha.getSize());
          sg::base::DataVector result_temp(alpha.getSize());
          sg::base::DataVector temp_two(alpha.getSize());

          #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp, result)
          {
            up(alpha, temp, this->algoDims[dim]);
            updown(temp, result, dim - 1, op_dim);
          }

          // Same from the other direction:
          #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp_two, result_temp)
          {
            updown(alpha, temp_two, dim - 1, op_dim);
            down(temp_two, result_temp, this->algoDims[dim]);
          }

          #pragma omp taskwait

          result.add(result_temp);
        } else {
          // Terminates dimension recursion
          sg::base::DataVector temp(alpha.getSize());

          #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, result)
          up(alpha, result, this->algoDims[dim]);

          #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp)
          down(alpha, temp, this->algoDims[dim]);

          #pragma omp taskwait

          result.add(temp);
        }

      }
    }

    void UpDownOneOpDim::specialOP(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim) {
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
          upOpDim(alpha, temp, this->algoDims[dim]);
          updown(temp, result, dim - 1, op_dim);
        }

        // Same from the other direction:
        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp_two, result_temp)
        {
          updown(alpha, temp_two, dim - 1, op_dim);
          downOpDim(temp_two, result_temp, this->algoDims[dim]);
        }

        #pragma omp taskwait

        result.add(result_temp);
      } else {
        // Terminates dimension recursion
        sg::base::DataVector temp(alpha.getSize());

        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, result)
        upOpDim(alpha, result, this->algoDims[dim]);

        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp)
        downOpDim(alpha, temp, this->algoDims[dim]);

        #pragma omp taskwait

        result.add(temp);
      }
    }

  }
}
