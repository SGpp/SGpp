/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/algorithm/StdUpDown.hpp"

namespace sg {
  namespace pde {

    StdUpDown::StdUpDown(sg::base::GridStorage* storage) : storage(storage), algoDims(storage->getAlgorithmicDimensions()), numAlgoDims_(storage->getAlgorithmicDimensions().size()) {
    }

    StdUpDown::~StdUpDown() {
    }

    void StdUpDown::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector beta(result.getSize());
      result.setAll(0.0);
      #pragma omp parallel
      {
        #pragma omp single nowait
        {
          this->updown(alpha, beta, this->numAlgoDims_ - 1);
        }
      }
      result.add(beta);
    }

    void StdUpDown::multParallelBuildingBlock(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector beta(result.getSize());
      result.setAll(0.0);

      this->updown(alpha, beta, this->numAlgoDims_ - 1);

      result.add(beta);
    }

    void StdUpDown::updown(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
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
          updown(temp, result, dim - 1);
        }


        #pragma omp task if(curNumAlgoDims - dim <= curMaxParallelDims) shared(alpha, temp_two, result_temp)
        {
          updown(alpha, temp_two, dim - 1);
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
}

