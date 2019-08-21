// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationDensitySamplingLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityConditional.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensitySampling1D.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/globaldef.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sgpp {
namespace datadriven {

void OperationDensitySamplingLinear::doSampling(base::DataVector* alpha, base::DataMatrix*& samples,
                                                size_t num_samples) {
  size_t num_dims = this->grid->getDimension();

  // output matrix
  samples = new base::DataMatrix(num_samples, num_dims);

  size_t size = num_samples / num_dims;

  if (size <= 0)
    throw base::operation_exception(
        "Error: # of dimensions greater than # of samples. Operation aborted!");

  size_t trunk = size;

  for (size_t dim_start = 0; dim_start < num_dims; dim_start++) {
    if (dim_start == num_dims - 1) size += num_samples % num_dims;

    // 1. marginalize to dim_start
    base::Grid* g1d = NULL;
    base::DataVector* a1d = NULL;
    op_factory::createOperationDensityMargTo1D(*this->grid)->margToDimX(alpha, g1d, a1d, dim_start);

    // 2. 1D sampling on dim_start
    base::DataVector* samples_start = NULL;
    unsigned int tseedp = static_cast<unsigned int>(static_cast<double>(time(NULL)) * 0.0001);
    op_factory::createOperationDensitySampling1D(*g1d)->doSampling1D(
        a1d, size, samples_start, &tseedp);
    delete g1d;
    delete a1d;

// 3. for every sample do...
#pragma omp parallel
    {
      base::DataVector* sampleVec = new base::DataVector(num_dims);
      unsigned int seedp;
      size_t samplesSize = samples_start->getSize();
#pragma omp critical
      {
#ifdef _WIN32
        double a = static_cast<double>(rand()) / RAND_MAX;
        int64_t b = static_cast<int64_t>(rand());
#else
        double a = static_cast<double>(rand_r(&tseedp)) / RAND_MAX;
        int64_t b = static_cast<int64_t>(rand_r(&tseedp));
#endif
#ifdef _OPENMP
        seedp =
            static_cast<unsigned int>(static_cast<double>(time(NULL)) * a +
                                      static_cast<double>((omp_get_thread_num() + 1) * 1000 * b));
#else
        seedp = static_cast<unsigned int>(static_cast<double>(time(NULL)) * a +
                                          static_cast<double>((1 + 1) * 1000 * b));
#endif
      }
#pragma omp for schedule(dynamic)

      for (size_t i = 0; i < samplesSize; i++) {
        sampleVec->set(dim_start, samples_start->get(i));
        doSampling_start_dimX(this->grid, alpha, dim_start, sampleVec, &seedp);

        for (size_t j = 0; j < num_dims; j++)
          samples->set(dim_start * trunk + i, j, sampleVec->get(j));
      }

      delete sampleVec;
    }
    delete samples_start;
  }

  return;
}

void OperationDensitySamplingLinear::doSampling(base::DataVector* alpha, base::DataMatrix*& samples,
                                                size_t num_samples, size_t dim_x) {
  size_t num_dims = this->grid->getDimension();

  if ((dim_x >= num_dims))
    throw base::operation_exception("Error: starting dimension out of range. Operation aborted!");

  // output matrix
  samples = new base::DataMatrix(num_samples, num_dims);

  // 1. marginalize to dim_x
  base::Grid* g1d = NULL;
  base::DataVector* a1d = NULL;
  op_factory::createOperationDensityMargTo1D(*this->grid)->margToDimX(alpha, g1d, a1d, dim_x);

  // 2. 1D sampling on dim_start
  base::DataVector* samples_start = NULL;
  unsigned int tseedp = static_cast<unsigned int>(static_cast<double>(time(NULL)) * 0.0001);
  op_factory::createOperationDensitySampling1D(*g1d)->doSampling1D(
      a1d, num_samples, samples_start, &tseedp);
  delete g1d;
  delete a1d;

// 3. for every sample do...
#pragma omp parallel
  {
    base::DataVector* sampleVec = new base::DataVector(num_dims);
    unsigned int seedp = 0;

#pragma omp critical
    {
#ifdef _WIN32
      double a = static_cast<double>(rand()) / RAND_MAX;
      int64_t b = static_cast<int64_t>(rand());
#else
      double a = static_cast<double>(rand_r(&tseedp)) / RAND_MAX;
      int64_t b = static_cast<int64_t>(rand_r(&tseedp));
#endif
#ifdef _OPENMP
      seedp =
          static_cast<unsigned int>(static_cast<double>(time(NULL)) * a +
                                    static_cast<double>((omp_get_thread_num() + 1) * 1000 * b));
#else
      seedp = static_cast<unsigned int>(static_cast<double>(time(NULL)) * a +
                                        static_cast<double>((1 + 1) * 1000 * b));
#endif
    }
#pragma omp for schedule(dynamic)

    for (size_t i = 0; i < num_samples; i++) {
      sampleVec->set(dim_x, samples_start->get(i));
      doSampling_start_dimX(this->grid, alpha, dim_x, sampleVec, &seedp);

      for (size_t j = 0; j < num_dims; j++) samples->set(i, j, sampleVec->get(j));
    }
  }
  delete samples_start;
  return;
}

void OperationDensitySamplingLinear::doSampling_start_dimX(base::Grid* g_in, base::DataVector* a_in,
                                                           size_t dim_start,
                                                           base::DataVector*& sampleVec,
                                                           unsigned int* seedp) {
  size_t dims = sampleVec->getSize();  // total dimensions

  if ((dims > 1) && (dim_start <= dims - 1)) {
    size_t curr_dim = dim_start;
    doSampling_in_next_dim(g_in, a_in, dim_start, sampleVec, curr_dim, seedp);
  } else if (dims == 1) {
    throw base::operation_exception("Error: # of dimensions = 1. No operation needed!");
  } else {
    throw base::operation_exception("Error: dimension out of range. Operation aborted!");
  }

  return;
}

void OperationDensitySamplingLinear::doSampling_in_next_dim(base::Grid* g_in,
                                                            base::DataVector* a_in, size_t dim_x,
                                                            base::DataVector*& sampleVec,
                                                            size_t& curr_dim, unsigned int* seedp) {
  size_t dims = sampleVec->getSize();  // total dimensions
  unsigned int op_dim =
      (curr_dim < dim_x) ? 0 : static_cast<unsigned int>(dim_x);  // actual dim to be operated on

  /* Step 1: do conditional in current dim */
  base::Grid* g_out = NULL;
  base::DataVector* a_out = new base::DataVector(1);
  op_factory::createOperationDensityConditional(*g_in)->doConditional(
      *a_in, g_out, *a_out, op_dim, sampleVec->get(curr_dim));

  // move on to next dim
  curr_dim = (curr_dim + 1) % dims;
  op_dim = (curr_dim < dim_x) ? 0 : static_cast<unsigned int>(dim_x);

  /* Step 2: draw a sample in next dim */
  base::DataVector* sample = NULL;

  if (g_out->getDimension() > 1) {
    // Marginalize to next dimension
    base::Grid* g1d = NULL;
    base::DataVector* a1d = NULL;
    op_factory::createOperationDensityMargTo1D(*g_out)->margToDimX(a_out, g1d, a1d, op_dim);

    // Draw a sample in next dimension
    op_factory::createOperationDensitySampling1D(*g1d)->doSampling1D(a1d, 1, sample, seedp);
    delete g1d;
    delete a1d;

  } else {
    // skip Marginalize, directly draw a sample in next dimension
    op_factory::createOperationDensitySampling1D(*g_out)->doSampling1D(a_out, 1, sample, seedp);
  }

  /* Step 3: copy sample to output */
  sampleVec->set(curr_dim, sample->get(0));
  delete sample;

  /* Step 4: sample in next dimension */
  if (g_out->getDimension() > 1)
    doSampling_in_next_dim(g_out, a_out, dim_x, sampleVec, curr_dim, seedp);

  delete g_out;
  delete a_out;

  return;
}
}  // namespace datadriven
}  // namespace sgpp
