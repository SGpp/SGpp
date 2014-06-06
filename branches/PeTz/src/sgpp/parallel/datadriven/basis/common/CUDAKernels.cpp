/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/common/operation/datadriven/CUDAKernels.hpp"
#include "tools/common/SGppStopwatch.hpp"

// including CUDA
#include <cuda.h>

#include <iostream>

float* gpu_grid_level_sp;
float* gpu_grid_index_sp;
float* gpu_dataset_sp;
float* gpu_alpha_sp;
float* gpu_datavec_sp;

__global__ void multTransSP_CUDA(float* ptrSource,
                                 float* ptrData,
                                 float* ptrLevel,
                                 float* ptrIndex,
                                 float* ptrResult,
                                 unsigned int sourceSize,
                                 unsigned int dims,
                                 unsigned int offset) {
  int globalIdx = (blockIdx.x * blockDim.x) + threadIdx.x + offset;
  int localIdx = threadIdx.x;
  float eval, index_calc, abs, last, localSupport, curSupport;
  float myResult = 0.0f;
  float locLevel[5];
  float locIndex[5];

  __shared__ float locData[320];
  __shared__ float locSource[64];

#pragma unroll 5

  for (unsigned int d = 0; d < dims; d++) {
    locLevel[d] = ptrLevel[(globalIdx * dims) + d];
    locIndex[d] = ptrIndex[(globalIdx * dims) + d];
  }

  // Iterate over all instances of the dataset
  for (unsigned int i = 0; i < sourceSize; i += 64) {
#pragma unroll 5

    for (unsigned int d = 0; d < dims; d++) {
      locData[(localIdx * dims) + d] = ptrData[((i + localIdx) * dims) + d];
    }

    locSource[localIdx] = ptrSource[i + localIdx];

    // Wait until all data is in shared memory
    __syncthreads();

    for (unsigned int k = 0; k < 64; k++) {
      curSupport = locSource[k];

#pragma unroll 5

      for (unsigned int d = 0; d < dims; d++) {
        eval = locLevel[d] * locData[(k * dims) + d];
        index_calc = eval - locIndex[d];
        abs = fabsf(index_calc);
        last = 1.0f - abs;
        localSupport = fmaxf(last, 0.0f);
        curSupport *= localSupport;
      }

      myResult += curSupport;
    }

    // Wait for calculation to finish (before new data can be loaded into shared memory)
    __syncthreads();
  }

  ptrResult[globalIdx] = myResult;
}

__global__ void multTransSP_CUDA_5d(float* ptrSource,
                                    float* ptrData,
                                    float* ptrLevel,
                                    float* ptrIndex,
                                    float* ptrResult,
                                    unsigned int sourceSize,
                                    unsigned int offset) {
  int globalIdx = (blockIdx.x * blockDim.x) + threadIdx.x + offset;
  int localIdx = threadIdx.x;
  float eval, index_calc, abs, last, localSupport, curSupport;
  float myResult = 0.0f;
  float locLevel_0;
  float locIndex_0;
  float locLevel_1;
  float locIndex_1;
  float locLevel_2;
  float locIndex_2;
  float locLevel_3;
  float locIndex_3;
  float locLevel_4;
  float locIndex_4;

  __shared__ float locData[320];
  __shared__ float locSource[64];

  locLevel_0 = ptrLevel[(globalIdx * 5) + 0];
  locIndex_0 = ptrIndex[(globalIdx * 5) + 0];
  locLevel_1 = ptrLevel[(globalIdx * 5) + 1];
  locIndex_1 = ptrIndex[(globalIdx * 5) + 1];
  locLevel_2 = ptrLevel[(globalIdx * 5) + 2];
  locIndex_2 = ptrIndex[(globalIdx * 5) + 2];
  locLevel_3 = ptrLevel[(globalIdx * 5) + 3];
  locIndex_3 = ptrIndex[(globalIdx * 5) + 3];
  locLevel_4 = ptrLevel[(globalIdx * 5) + 4];
  locIndex_4 = ptrIndex[(globalIdx * 5) + 4];

  // Iterate over all instances of the dataset
  for (unsigned int i = 0; i < sourceSize; i += 64) {
    locData[(localIdx * 5) + 0] = ptrData[((i + localIdx) * 5) + 0];
    locData[(localIdx * 5) + 1] = ptrData[((i + localIdx) * 5) + 1];
    locData[(localIdx * 5) + 2] = ptrData[((i + localIdx) * 5) + 2];
    locData[(localIdx * 5) + 3] = ptrData[((i + localIdx) * 5) + 3];
    locData[(localIdx * 5) + 4] = ptrData[((i + localIdx) * 5) + 4];

    locSource[localIdx] = ptrSource[i + localIdx];

    // Wait until all data is in shared memory
    __syncthreads();

    for (unsigned int k = 0; k < 64; k++) {
      curSupport = locSource[k];

      eval = locLevel_0 * locData[(k * 5) + 0];
      index_calc = eval - locIndex_0;
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel_1 * locData[(k * 5) + 1];
      index_calc = eval - locIndex_1;
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel_2 * locData[(k * 5) + 2];
      index_calc = eval - locIndex_2;
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel_3 * locData[(k * 5) + 3];
      index_calc = eval - locIndex_3;
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel_4 * locData[(k * 5) + 4];
      index_calc = eval - locIndex_4;
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      myResult += curSupport;
    }

    // Wait for calculation to finish (before new data can be loaded into shared memory)
    __syncthreads();
  }

  ptrResult[globalIdx] = myResult;
}

__global__ void multSP_CUDA(float* ptrAlpha,
                            float* ptrData,
                            float* ptrLevel,
                            float* ptrIndex,
                            float* ptrResult,
                            unsigned int fastStorageSize,
                            unsigned int storageSize,
                            unsigned int dims,
                            unsigned int offset) {
  int globalIdx = (blockIdx.x * blockDim.x) + threadIdx.x + offset;
  int localIdx = threadIdx.x;
  float eval, index_calc, abs, last, localSupport, curSupport;
  float myResult = 0.0f;
  float locData[5];

  __shared__ float locLevel[320];
  __shared__ float locIndex[320];
  __shared__ float locAlpha[64];

#pragma unroll 5

  for (unsigned int d = 0; d < dims; d++) {
    locData[d] = ptrData[(globalIdx * dims) + d];
  }

  // Iterate over all grid points (fast ones, with cache)
  for (unsigned int j = 0; j < fastStorageSize; j += 64) {
#pragma unroll 5

    for (unsigned int d = 0; d < dims; d++) {
      locLevel[(localIdx * dims) + d] = ptrLevel[((j + localIdx) * dims) + d];
      locIndex[(localIdx * dims) + d] = ptrIndex[((j + localIdx) * dims) + d];
    }

    locAlpha[localIdx] = ptrAlpha[j + localIdx];

    // Wait until all needed data is in shared memory
    __syncthreads();

    for (unsigned int k = 0; k < 64; k++) {
      curSupport = locAlpha[k];

#pragma unroll 5

      for (unsigned int d = 0; d < dims; d++) {
        eval = locLevel[(k * dims) + d] * locData[d];
        index_calc = eval - locIndex[(k * dims) + d];
        abs = fabsf(index_calc);
        last = 1.0f - abs;
        localSupport = fmaxf(last, 0.0f);
        curSupport *= localSupport;
      }

      myResult += curSupport;
    }

    // wait until calculations have finished (avoid overwriting shared memory)
    __syncthreads();
  }

  // Iterate over all grid points (slow ones, without cache)
  for (unsigned int m = fastStorageSize; m < storageSize; m++) {
    curSupport = ptrAlpha[m];

    for (unsigned int d = 0; d < dims; d++) {
      eval = ptrLevel[(m * dims) + d] * locData[d];
      index_calc = eval - ptrIndex[(m * dims) + d];
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;
    }

    myResult += curSupport;
  }

  ptrResult[globalIdx] = myResult;
}

__global__ void multSP_CUDA_5d(float* ptrAlpha,
                               float* ptrData,
                               float* ptrLevel,
                               float* ptrIndex,
                               float* ptrResult,
                               unsigned int fastStorageSize,
                               unsigned int storageSize,
                               unsigned int offset) {
  int globalIdx = (blockIdx.x * blockDim.x) + threadIdx.x + offset;
  int localIdx = threadIdx.x;
  float eval, index_calc, abs, last, localSupport, curSupport;
  float myResult = 0.0f;
  float locData_0;
  float locData_1;
  float locData_2;
  float locData_3;
  float locData_4;

  __shared__ float locLevel[320];
  __shared__ float locIndex[320];
  __shared__ float locAlpha[64];

  locData_0 = ptrData[(globalIdx * 5) + 0];
  locData_1 = ptrData[(globalIdx * 5) + 1];
  locData_2 = ptrData[(globalIdx * 5) + 2];
  locData_3 = ptrData[(globalIdx * 5) + 3];
  locData_4 = ptrData[(globalIdx * 5) + 4];

  // Iterate over all grid points (fast ones, with cache)
  for (unsigned int j = 0; j < fastStorageSize; j += 64) {
    locLevel[(localIdx * 5) + 0] = ptrLevel[((j + localIdx) * 5) + 0];
    locIndex[(localIdx * 5) + 0] = ptrIndex[((j + localIdx) * 5) + 0];
    locLevel[(localIdx * 5) + 1] = ptrLevel[((j + localIdx) * 5) + 1];
    locIndex[(localIdx * 5) + 1] = ptrIndex[((j + localIdx) * 5) + 1];
    locLevel[(localIdx * 5) + 2] = ptrLevel[((j + localIdx) * 5) + 2];
    locIndex[(localIdx * 5) + 2] = ptrIndex[((j + localIdx) * 5) + 2];
    locLevel[(localIdx * 5) + 3] = ptrLevel[((j + localIdx) * 5) + 3];
    locIndex[(localIdx * 5) + 3] = ptrIndex[((j + localIdx) * 5) + 3];
    locLevel[(localIdx * 5) + 4] = ptrLevel[((j + localIdx) * 5) + 4];
    locIndex[(localIdx * 5) + 4] = ptrIndex[((j + localIdx) * 5) + 4];

    locAlpha[localIdx] = ptrAlpha[j + localIdx];

    // Wait until all needed data is in shared memory
    __syncthreads();

    for (unsigned int k = 0; k < 64; k++) {
      curSupport = locAlpha[k];

      eval = locLevel[(k * 5) + 0] * locData_0;
      index_calc = eval - locIndex[(k * 5) + 0];
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel[(k * 5) + 1] * locData_1;
      index_calc = eval - locIndex[(k * 5) + 1];
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel[(k * 5) + 2] * locData_2;
      index_calc = eval - locIndex[(k * 5) + 2];
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel[(k * 5) + 3] * locData_3;
      index_calc = eval - locIndex[(k * 5) + 3];
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel[(k * 5) + 4] * locData_4;
      index_calc = eval - locIndex[(k * 5) + 4];
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      myResult += curSupport;
    }

    // wait until calculations have finished (avoid overwriting shared memory)
    __syncthreads();
  }

  // Iterate over all grid points (slow ones, without cache)
  for (unsigned int m = fastStorageSize; m < storageSize; m++) {
    curSupport = ptrAlpha[m];

    eval = ptrLevel[(m * 5) + 0] * locData_0;
    index_calc = eval - ptrIndex[(m * 5) + 0];
    abs = fabsf(index_calc);
    last = 1.0f - abs;
    localSupport = fmaxf(last, 0.0f);
    curSupport *= localSupport;

    eval = ptrLevel[(m * 5) + 1] * locData_1;
    index_calc = eval - ptrIndex[(m * 5) + 1];
    abs = fabsf(index_calc);
    last = 1.0f - abs;
    localSupport = fmaxf(last, 0.0f);
    curSupport *= localSupport;

    eval = ptrLevel[(m * 5) + 2] * locData_2;
    index_calc = eval - ptrIndex[(m * 5) + 2];
    abs = fabsf(index_calc);
    last = 1.0f - abs;
    localSupport = fmaxf(last, 0.0f);
    curSupport *= localSupport;

    eval = ptrLevel[(m * 5) + 3] * locData_3;
    index_calc = eval - ptrIndex[(m * 5) + 3];
    abs = fabsf(index_calc);
    last = 1.0f - abs;
    localSupport = fmaxf(last, 0.0f);
    curSupport *= localSupport;

    eval = ptrLevel[(m * 5) + 4] * locData_4;
    index_calc = eval - ptrIndex[(m * 5) + 4];
    abs = fabsf(index_calc);
    last = 1.0f - abs;
    localSupport = fmaxf(last, 0.0f);
    curSupport *= localSupport;

    myResult += curSupport;
  }

  ptrResult[globalIdx] = myResult;
}

double multTransSPCUDA(float* ptrSource, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims) {
  // allocate memory on GPU
  cudaMalloc((void**) &gpu_alpha_sp, storageSize * sizeof(float));
  cudaMalloc((void**) &gpu_datavec_sp, sourceSize * sizeof(float));

  // copy coefficients to GPU
  cudaMemcpy(gpu_datavec_sp, ptrSource, sourceSize * sizeof(float), cudaMemcpyHostToDevice);
  int myStorageSize = ((int)storageSize) / 64;

  // Measure time
  sg::base::SGppStopwatch* mytimer = new sg::base::SGppStopwatch();
  mytimer->start();

  if (dims == 5) {
    multTransSP_CUDA_5d <<< myStorageSize, 64>>>(gpu_datavec_sp, gpu_dataset_sp, gpu_grid_level_sp, gpu_grid_index_sp, gpu_alpha_sp, (unsigned int)sourceSize, 0);
  } else {
    multTransSP_CUDA <<< myStorageSize, 64>>>(gpu_datavec_sp, gpu_dataset_sp, gpu_grid_level_sp, gpu_grid_index_sp, gpu_alpha_sp, (unsigned int)sourceSize, (unsigned int)dims, 0);
  }

  // copy results back to host
  cudaMemcpy(ptrGlobalResult, gpu_alpha_sp, storageSize * sizeof(float), cudaMemcpyDeviceToHost);

  double time = mytimer->stop();
  delete mytimer;

  // free data on GPU
  cudaFree(gpu_datavec_sp);
  cudaFree(gpu_alpha_sp);

  return time;
}

double multSPCUDA(float* ptrAlpha, float* ptrResult, size_t result_size, size_t storageSize, size_t dims) {
  // allocate memory on GPU
  cudaMalloc((void**) &gpu_alpha_sp, storageSize * sizeof(float));
  cudaMalloc((void**) &gpu_datavec_sp, result_size * sizeof(float));

  // copy coefficients to GPU
  cudaMemcpy(gpu_alpha_sp, ptrAlpha, storageSize * sizeof(float), cudaMemcpyHostToDevice);

  unsigned int tmp = ((unsigned int)storageSize) / CUDA_BLOCK_SIZE_GPU;
  unsigned int fastStorageSize = tmp * CUDA_BLOCK_SIZE_GPU;

  // Measure time
  sg::base::SGppStopwatch* mytimer = new sg::base::SGppStopwatch();
  mytimer->start();

  if (dims == 5) {
    multSP_CUDA_5d <<< result_size / 64, 64 >>> (gpu_alpha_sp, gpu_dataset_sp, gpu_grid_level_sp, gpu_grid_index_sp, gpu_datavec_sp, fastStorageSize, (unsigned int)storageSize, 0);
  } else {
    multSP_CUDA <<< result_size / 64, 64 >>> (gpu_alpha_sp, gpu_dataset_sp, gpu_grid_level_sp, gpu_grid_index_sp, gpu_datavec_sp, fastStorageSize, (unsigned int)storageSize, (unsigned int)dims, 0);
  }

  // copy results back to host
  cudaMemcpy(ptrResult, gpu_datavec_sp, result_size * sizeof(float), cudaMemcpyDeviceToHost);

  double time = mytimer->stop();
  delete mytimer;

  // free data on GPU
  cudaFree(gpu_datavec_sp);
  cudaFree(gpu_alpha_sp);

  return time;
}

void uploadGridSPCUDA(float* ptrLevel, float* ptrIndex, size_t storageSize, size_t dims) {
  size_t mem_size = storageSize * dims * sizeof(float);

  // allocate memory on GPU
  cudaMalloc((void**) &gpu_grid_level_sp, mem_size);
  cudaMalloc((void**) &gpu_grid_index_sp, mem_size);

  // copy host memory to device
  cudaMemcpy(gpu_grid_level_sp, ptrLevel, mem_size, cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_grid_index_sp, ptrIndex, mem_size, cudaMemcpyHostToDevice);
}

void uploadDataSPCUDA(float* ptrData, size_t dataSize, size_t dims) {
  size_t mem_size = dataSize * dims * sizeof(float);

  // allocate memory on GPU
  cudaMalloc((void**) &gpu_dataset_sp, mem_size);

  // copy host memory to device
  cudaMemcpy(gpu_dataset_sp, ptrData, mem_size, cudaMemcpyHostToDevice);
}

void deleteGridSPCUDA() {
  cudaFree(gpu_grid_level_sp);
  cudaFree(gpu_grid_index_sp);
}

void deleteDataSPCUDA() {
  cudaFree(gpu_dataset_sp);
}
