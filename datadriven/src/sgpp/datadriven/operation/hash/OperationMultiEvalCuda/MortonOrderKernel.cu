#include "basicCuda.hpp"
#include "cudaHelper.hpp"
#include "MortonOrderKernel.hpp"
#include "kernels.cuh"

#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/sort.h>
#include <thrust/copy.h>

/// Wrapper for the kernel call for Morton order curve index computation
void zorder (double* pos, size_t* perm, uint32_t m, uint32_t DIM) {
  initCudaConstants(DIM);
  CudaCheckError();
  thrust::device_ptr<size_t> perm_p(perm);
  thrust::device_vector<uint64_t> idx(m);
  // Compute index in the Morton order curve
  gpu_zindex<<<m/CUDA_BLOCKSIZE,CUDA_BLOCKSIZE>>>(pos, thrust::raw_pointer_cast(idx.data()));
  CudaCheckError();
  // Sort according to this index
  thrust::stable_sort_by_key(idx.begin(), idx.end(), perm_p);
}

