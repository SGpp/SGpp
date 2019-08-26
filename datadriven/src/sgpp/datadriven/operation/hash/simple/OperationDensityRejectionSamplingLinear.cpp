// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationDensityRejectionSamplingLinear.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {
void OperationDensityRejectionSamplingLinear::doSampling(base::DataVector* alpha,
                                                         base::DataMatrix*& samples,
                                                         size_t num_samples, size_t trial_max) {
  size_t num_dims = this->grid->getDimension();
  samples = new base::DataMatrix(num_samples, num_dims);  // output samples

  size_t SEARCH_MAX = 100000;  // find the approximated maximum of function with 100000 points
  double maxValue = 0;        // the approximated maximum value of function

  // search for (approx.) maximum of function
  base::DataMatrix* tmp = new base::DataMatrix(SEARCH_MAX, num_dims);
  base::DataVector* tmpEval = new base::DataVector(SEARCH_MAX);

#pragma omp parallel
  {
#ifndef _WIN32
#ifdef _OPENMP
    unsigned int seedp = static_cast<unsigned int>(static_cast<double>(time(nullptr)) *
                                                   static_cast<double>(omp_get_thread_num() + 1));
#else
    unsigned int seedp =
        static_cast<unsigned int>(static_cast<double>(time(nullptr)) * static_cast<double>(1 + 1));
#endif
#endif
#pragma omp for

    for (size_t i = 0; i < SEARCH_MAX; i++) {
      for (size_t j = 0; j < num_dims; j++)
#ifdef _WIN32
        tmp->set(i, j, static_cast<double>(rand()) / RAND_MAX);

#else
        tmp->set(i, j, static_cast<double>(rand_r(&seedp)) / RAND_MAX);
#endif
    }
  }

  op_factory::createOperationMultipleEval(*grid, *(tmp))->mult(*alpha, *tmpEval);
  maxValue = tmpEval->max();
  delete tmp;
  tmp = nullptr;
  delete tmpEval;
  tmpEval = nullptr;

#pragma omp parallel
  {
#ifndef _WIN32
#ifdef _OPENMP
    unsigned int seedp = static_cast<unsigned int>(time(nullptr)) * (omp_get_thread_num() + 1);
#else
    unsigned int seedp = static_cast<unsigned int>(time(nullptr)) * (1 + 1);
#endif
#endif
    base::DataVector p(num_dims);
    double fhat = 0.0;
    std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(*grid));

#pragma omp for schedule(dynamic)

    for (size_t i = 0; i < num_samples; i++) {  // for every sample
      // find the appropriate sample within a # of trials
      size_t j = 0;

      for (; j < trial_max; j++) {
        // pick a random data point "p"
        for (size_t d = 0; d < num_dims; d++)
#ifdef _WIN32
          p[d] = static_cast<double>(rand()) / RAND_MAX;

#else
          p[d] = static_cast<double>(rand_r(&seedp)) / RAND_MAX;
#endif

        // evaluate at this point "p"
        fhat = opEval->eval(*alpha, p);

#ifdef _WIN32

        if ((static_cast<double>(rand()) / RAND_MAX * maxValue < fhat) &&
            (fhat > maxValue * 0.01)) {
#else

        if ((static_cast<double>(rand_r(&seedp)) / RAND_MAX * maxValue < fhat) &&
            (fhat > maxValue * 0.01)) {
#endif
          samples->setRow(i, p);
          break;
        }
      }

      if (j == trial_max)
        throw base::operation_exception("Error: maximum # of trials reached. Operation aborted!");
    }
  }

  return;
}  // end of doSampling()
}  // namespace datadriven
}  // namespace sgpp
