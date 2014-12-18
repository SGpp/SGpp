/* *****************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                          *
* This file is part of the SG++ project. For conditions of distribution and    *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp           *
***************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "misc/pde/algorithm/UpDownOneOpDimEnhanced.hpp"

namespace sg {
  namespace pde {

    UpDownOneOpDimEnhanced::UpDownOneOpDimEnhanced(sg::base::GridStorage* storage, sg::base::DataVector& coef) : storage(storage), coefs(&coef), algoDims(storage->getAlgorithmicDimensions()), numAlgoDims_(storage->getAlgorithmicDimensions().size()) {
    }

    UpDownOneOpDimEnhanced::UpDownOneOpDimEnhanced(sg::base::GridStorage* storage) : storage(storage), coefs(NULL), algoDims(storage->getAlgorithmicDimensions()), numAlgoDims_(storage->getAlgorithmicDimensions().size()) {
    }

    UpDownOneOpDimEnhanced::~UpDownOneOpDimEnhanced() {
    }

    void UpDownOneOpDimEnhanced::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataMatrix beta(result.getSize(), this->numAlgoDims_);
      sg::base::DataMatrix maAlpha(alpha.getSize(), this->numAlgoDims_);

      result.setAll(0.0);
      maAlpha.expand(alpha);

      #pragma omp parallel
      {
        #pragma omp single nowait
        {
          this->updown(maAlpha, beta, this->numAlgoDims_ - 1);
        }
      }

      if (coefs == NULL) {
        beta.addReduce(result);
      } else {
        beta.addReduce(result, *coefs, 0);
      }
    }

    void UpDownOneOpDimEnhanced::multParallelBuildingBlock(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataMatrix beta(result.getSize(), this->numAlgoDims_);
      sg::base::DataMatrix maAlpha(alpha.getSize(), this->numAlgoDims_);

      result.setAll(0.0);
      maAlpha.expand(alpha);

      this->updown(maAlpha, beta, this->numAlgoDims_ - 1);

      if (coefs == NULL) {
        beta.addReduce(result);
      } else {
        beta.addReduce(result, *coefs, 0);
      }
    }

    void UpDownOneOpDimEnhanced::updown(sg::base::DataMatrix& alpha, sg::base::DataMatrix& result, size_t dim) {
      size_t curNumAlgoDims = this->numAlgoDims_;
      size_t curMaxParallelDims = this->maxParallelDims_;

      //Unidirectional scheme
      if (dim > 0) {
        // Reordering ups and downs
        sg::base::DataMatrix temp(alpha.getNrows(), this->numAlgoDims_);
        sg::base::DataMatrix result_temp(alpha.getNrows(), this->numAlgoDims_);
        sg::base::DataMatrix temp_two(alpha.getNrows(), this->numAlgoDims_);

        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp, result)
        {
          up(alpha, temp, dim);
          updown(temp, result, dim - 1);
        }

        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp_two, result_temp)
        {
          updown(alpha, temp_two, dim - 1);
          down(temp_two, result_temp, dim);
        }

        #pragma omp taskwait

        result.add(result_temp);
      } else {
        // Terminates dimension recursion
        sg::base::DataMatrix temp(alpha.getNrows(), this->numAlgoDims_);

        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, result)
        up(alpha, result, dim);

        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp)
        down(alpha, temp, dim);

        #pragma omp taskwait

        result.add(temp);
      }
    }

  }
}
