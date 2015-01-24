/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/algorithm/UpDownTwoOpDims.hpp"

namespace sg {
  namespace pde {

    UpDownTwoOpDims::UpDownTwoOpDims(sg::base::GridStorage* storage, sg::base::DataMatrix& coef) : storage(storage), coefs(&coef), algoDims(storage->getAlgorithmicDimensions()), numAlgoDims_(storage->getAlgorithmicDimensions().size()) {
    }

    UpDownTwoOpDims::UpDownTwoOpDims(sg::base::GridStorage* storage) : storage(storage), coefs(NULL), algoDims(storage->getAlgorithmicDimensions()), numAlgoDims_(storage->getAlgorithmicDimensions().size()) {
    }

    UpDownTwoOpDims::~UpDownTwoOpDims() {
    }

    void UpDownTwoOpDims::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      result.setAll(0.0);

      #pragma omp parallel
      {
        #pragma omp single nowait
        {
          for (size_t i = 0; i < this->numAlgoDims_; i++) {
            for (size_t j = 0; j < this->numAlgoDims_; j++) {
              // use the operator's symmetry
              if ( j <= i) {
                #pragma omp task firstprivate(i, j) shared(alpha, result)
                {
                  sg::base::DataVector beta(result.getSize());

                  if (this->coefs != NULL) {
                    if (this->coefs->get(i, j) != 0.0) {
                      this->updown(alpha, beta, this->numAlgoDims_ - 1, i, j);

                      #pragma omp critical
                      {
                        result.axpy(this->coefs->get(i, j), beta);
                      }
                    }
                  } else
                  {
                    this->updown(alpha, beta, this->numAlgoDims_ - 1, i, j);

                    #pragma omp critical
                    {
                      result.add(beta);
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

    void UpDownTwoOpDims::multParallelBuildingBlock(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t operationDimOne, size_t operationDimTwo) {
      result.setAll(0.0);
      sg::base::DataVector beta(result.getSize());

      // use the operator's symmetry
      if ( operationDimTwo <= operationDimOne) {
        if (this->coefs != NULL) {
          if (this->coefs->get(operationDimOne, operationDimTwo) != 0.0) {
            this->updown(alpha, beta, this->numAlgoDims_ - 1, operationDimOne, operationDimTwo);
            result.axpy(this->coefs->get(operationDimOne, operationDimTwo), beta);
          }
        } else {
          this->updown(alpha, beta, this->numAlgoDims_ - 1, operationDimOne, operationDimTwo);
          result.add(beta);
        }
      }
    }

    void UpDownTwoOpDims::updown(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two) {
      if ((dim == op_dim_one) && (dim == op_dim_two)) {
        specialOpOneAndOpTwo(alpha, result, dim, op_dim_one, op_dim_two);
      } else if ((dim == op_dim_one || dim == op_dim_two) && (op_dim_one != op_dim_two)) {
        if (dim == op_dim_one) {
          specialOpOne(alpha, result, dim, op_dim_one, op_dim_two);
        }

        if (dim == op_dim_two) {
          specialOpTwo(alpha, result, dim, op_dim_one, op_dim_two);
        }
      } else {
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
            up(alpha, temp, this->algoDims[dim]);
            updown(temp, result, dim - 1, op_dim_one, op_dim_two);
          }

          // Same from the other direction:
          #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp_two, result_temp)
          {
            updown(alpha, temp_two, dim - 1, op_dim_one, op_dim_two);
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

    void UpDownTwoOpDims::specialOpOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two) {
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
          upOpDimOne(alpha, temp, this->algoDims[dim]);
          updown(temp, result, dim - 1, op_dim_one, op_dim_two);
        }

        // Same from the other direction:
        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp_two, result_temp)
        {
          updown(alpha, temp_two, dim - 1, op_dim_one, op_dim_two);
          downOpDimOne(temp_two, result_temp, this->algoDims[dim]);
        }

        #pragma omp taskwait

        result.add(result_temp);
      } else {
        // Terminates dimension recursion
        sg::base::DataVector temp(alpha.getSize());

        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, result)
        upOpDimOne(alpha, result, this->algoDims[dim]);

        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp)
        downOpDimOne(alpha, temp, this->algoDims[dim]);

        #pragma omp taskwait

        result.add(temp);
      }
    }

    void UpDownTwoOpDims::specialOpTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two) {
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
          upOpDimTwo(alpha, temp, this->algoDims[dim]);
          updown(temp, result, dim - 1, op_dim_one, op_dim_two);
        }

        // Same from the other direction:
        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp_two, result_temp)
        {
          updown(alpha, temp_two, dim - 1, op_dim_one, op_dim_two);
          downOpDimTwo(temp_two, result_temp, this->algoDims[dim]);
        }

        #pragma omp taskwait

        result.add(result_temp);
      } else {
        // Terminates dimension recursion
        sg::base::DataVector temp(alpha.getSize());

        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, result)
        upOpDimTwo(alpha, result, this->algoDims[dim]);

        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp)
        downOpDimTwo(alpha, temp, this->algoDims[dim]);

        #pragma omp taskwait

        result.add(temp);
      }
    }

    void UpDownTwoOpDims::specialOpOneAndOpTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two) {
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
          upOpDimOneAndOpDimTwo(alpha, temp, this->algoDims[dim]);
          updown(temp, result, dim - 1, op_dim_one, op_dim_two);
        }

        // Same from the other direction:
        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp_two, result_temp)
        {
          updown(alpha, temp_two, dim - 1, op_dim_one, op_dim_two);
          downOpDimOneAndOpDimTwo(temp_two, result_temp, this->algoDims[dim]);
        }

        #pragma omp taskwait

        result.add(result_temp);
      } else {
        // Terminates dimension recursion
        sg::base::DataVector temp(alpha.getSize());

        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, result)
        upOpDimOneAndOpDimTwo(alpha, result, this->algoDims[dim]);

        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp)
        downOpDimOneAndOpDimTwo(alpha, temp, this->algoDims[dim]);

        #pragma omp taskwait

        result.add(temp);
      }
    }

  }
}
