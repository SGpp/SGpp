#include "basicCuda.hpp"
#include "cudaHelper.hpp"
#include "kernels.cuh"
#include "MultiEvalKernel.hpp"
#include "consts.hpp"

/// Wrapper for kernel call of the stream boundary limitation
void streamboundCuda (double* pos, gridnode_t* node, limit_t* limit, uint32_t M, uint32_t N) {
  uint64_t *idx_p;
	cudaMalloc((void**)&idx_p,sizeof(uint64_t)*M);
  cudaDeviceSynchronize();
  CudaCheckError();
  // Compute index of the Morton order curve
  gpu_zindex<<<M/CUDA_BLOCKSIZE,CUDA_BLOCKSIZE>>>(pos, idx_p);
 	cudaDeviceSynchronize();
  CudaCheckError();
  // Compute lower and upper boundaries
  gpu_zbound<<<N/CUDA_BLOCKSIZE,CUDA_BLOCKSIZE>>>(idx_p, node, limit, M);
 	cudaDeviceSynchronize();
  CudaCheckError();
  cudaFree(idx_p);
}

/// Wrapper for kernel call of subspace limits for standard eval
void preprocessCuda (gridnode_t* node, double* pos, uint32_t* limit, uint32_t maxlevel, uint32_t M,
    uint32_t DIM) {
	initCudaConstants(DIM);
  gpu_preprocess<<<(M/CUDA_BLOCKSIZE),CUDA_BLOCKSIZE>>>(limit, node, pos, maxlevel);
  cudaDeviceSynchronize();
}

/// Wrapper for kernel call of standard evaluation
void evalCuda (double* res, double *a, gridnode_t* node, double* pos, uint32_t M,
  uint32_t maxlevel, uint32_t* limit, uint32_t subcnt, uint32_t* subs) {
  gpu_eval<<<(M/CUDA_BLOCKSIZE),CUDA_BLOCKSIZE>>>(res, node, a, pos, limit, subcnt, subs);
  cudaDeviceSynchronize();
}

/// Wrapper for kernel call of transposed eval with additional FMA
void transposedCuda (double* a, gridnode_t* node, double* pos, double* y,
    limit_t* limit, double* b, double c, uint32_t M, uint32_t N) {
  gpu_transevel<<<N,CUDA_BLOCKSIZE>>>(a, node, pos, y, limit, b, c, 1.0/double(M));
  cudaDeviceSynchronize();
}

/// Wrapper for kernel call of transposed eval
void transposedCuda (double* a, gridnode_t* node, double* pos, double* y,
    limit_t* limit, uint32_t M, uint32_t N) {
  gpu_transevel<<<N,CUDA_BLOCKSIZE>>>(a, node, pos, y, limit);
  cudaDeviceSynchronize();
}
