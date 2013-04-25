/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author A. Mo-Hellenbrand

#include "datadriven/basis/linear/noboundary/operation/OperationDensityRejectionSamplingLinear.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "datadriven/operation/DatadrivenOpFactory.hpp"
#include "base/exception/operation_exception.hpp"
#include <omp.h>

namespace sg {
  namespace datadriven {
    void OperationDensityRejectionSamplingLinear::doSampling(base::DataVector* alpha, base::DataMatrix*& samples, size_t num_samples, size_t trial_max) {

      size_t num_dims = this->grid->getStorage()->dim();
      samples = new base::DataMatrix(num_samples, num_dims); //output samples

      size_t SEARCH_MAX = 100000; //find the approximated maximum of function with 100000 points
      double maxValue = 0; //the approximated maximum value of function

      //search for (approx.) maximum of function
      base::DataMatrix* tmp = new base::DataMatrix(SEARCH_MAX, num_dims);
      base::DataVector* tmpEval = new base::DataVector(SEARCH_MAX);

      for (size_t i = 0; i < SEARCH_MAX; i++) {
        for (size_t j = 0; j < num_dims; j++)
          tmp->set(i, j, (double)rand() / RAND_MAX);
      }

      base::OperationMultipleEval* opMultEval = op_factory::createOperationMultipleEval(*grid, tmp);
      opMultEval->mult(*alpha, *tmpEval);
      maxValue = tmpEval->max();
      delete tmp;
      tmp = NULL;
      delete tmpEval;
      tmpEval = NULL;

      double fhat = 0.0;
      size_t j;
      base::DataVector* p;
      base::OperationEval* opEval;

      #pragma omp parallel for private(p, j, opEval) schedule(dynamic)

      for (size_t i = 0; i < num_samples; i++) { //for every sample

        p = new base::DataVector(num_dims);
        opEval = op_factory::createOperationEval(*grid);

        //find the appropriate sample within a # of trial
        for (j = 0; j < trial_max; j++) {

          // pick a random data point "p"
          for (size_t d = 0; d < num_dims; d++)
            p->set(d, (double)rand() / RAND_MAX);

          // evaluate at this point "p"
          fhat = opEval->eval(*alpha, *p);

          if (((double)rand() / RAND_MAX * maxValue < fhat) && (fhat > maxValue * 0.050)) {
            samples->setRow(i, *p);
            break;
          }
        }

        if (j == trial_max)
          throw base::operation_exception("Error: maximum # of trials reached. Operation aborted!");

        delete p;
        delete opEval;
      }

      return;
    } //end of doSampling()

  }
}
