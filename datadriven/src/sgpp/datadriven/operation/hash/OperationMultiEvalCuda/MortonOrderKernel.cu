// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <thrust/copy.h>
#include <thrust/device_malloc.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

#include "MortonOrderKernel.hpp"
#include "basicCuda.hpp"
#include "cudaHelper.hpp"
#include "kernels.cuh"

///@cond DOXY_IGNORE // NOLINT()
namespace sgpp {
namespace datadriven {
namespace OpMultiEvalCudaDetail {

/// Wrapper for the kernel call for Morton order curve index computation
void zorder(double* pos, size_t* perm, size_t m, size_t DIM) {
  initCudaConstants(static_cast<uint32_t>(DIM));
  CudaCheckError();
  thrust::device_ptr<size_t> perm_p(perm);
  thrust::device_vector<uint64_t> idx(m);
  // Compute index in the Morton order curve
  gpu_zindex<<<m / CUDA_BLOCKSIZE, CUDA_BLOCKSIZE>>>(pos, thrust::raw_pointer_cast(idx.data()));
  CudaCheckError();
  // Sort according to this index
  thrust::stable_sort_by_key(idx.begin(), idx.end(), perm_p);
}

}  // namespace OpMultiEvalCudaDetail
}  // namespace datadriven
}  // namespace sgpp
///@endcond // NOLINT()
