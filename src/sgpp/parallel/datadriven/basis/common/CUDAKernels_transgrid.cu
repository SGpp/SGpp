/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/tools/SGppStopwatch.hpp"

#include "parallel/datadriven/basis/common/CUDAKernels.hpp"

// including CUDA
#include <cuda.h>

#include <iostream>

float* gpu_grid_level_sp;
float* gpu_grid_index_sp;
float* gpu_dataset_sp;
float* gpu_alpha_sp;
float* gpu_datavec_sp;
float* host_alpha_sp;
float* host_datavec_sp;
float* host_grid_level_sp;
float* host_grid_index_sp;
size_t gpu_full_storageSize;


__global__ void multTransSP_CUDA(float* ptrSource,
                                 float* ptrData,
                                 float* ptrLevel,
                                 float* ptrIndex,
                                 float* ptrResult,
                                 unsigned int sourceSize,
                                 unsigned int storageSize,
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
    locLevel[d] = ptrLevel[globalIdx + (storageSize * d)];
    locIndex[d] = ptrIndex[globalIdx + (storageSize * d)];
  }

  // Iterate over all instances of the dataset
  for (unsigned int i = 0; i < sourceSize; i += 64) {
#pragma unroll 5

    for (unsigned int d = 0; d < dims; d++) {
      locData[(d*64)+localIdx] = ptrData[(d*sourceSize)+(i+localIdx)];
    }

    locSource[localIdx] = ptrSource[i + localIdx];

    // Wait until all data is in shared memory
    __syncthreads();

    for (unsigned int k = 0; k < 64; k++) {
      curSupport = locSource[k];

#pragma unroll 5

      for (unsigned int d = 0; d < dims; d++) {
        eval = locLevel[d] * locData[(d*64) + k];
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
                                    unsigned int storageSize,
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

  locLevel_0 = ptrLevel[globalIdx + (storageSize * 0)];
  locIndex_0 = ptrIndex[globalIdx + (storageSize * 0)];
  locLevel_1 = ptrLevel[globalIdx + (storageSize * 1)];
  locIndex_1 = ptrIndex[globalIdx + (storageSize * 1)];
  locLevel_2 = ptrLevel[globalIdx + (storageSize * 2)];
  locIndex_2 = ptrIndex[globalIdx + (storageSize * 2)];
  locLevel_3 = ptrLevel[globalIdx + (storageSize * 3)];
  locIndex_3 = ptrIndex[globalIdx + (storageSize * 3)];
  locLevel_4 = ptrLevel[globalIdx + (storageSize * 4)];
  locIndex_4 = ptrIndex[globalIdx + (storageSize * 4)];

  // Iterate over all instances of the dataset
  for (unsigned int i = 0; i < sourceSize; i += 64) {
    locData[(0*64)+localIdx] = ptrData[(i + localIdx) + (0*sourceSize)];
    locData[(1*64)+localIdx] = ptrData[(i + localIdx) + (1*sourceSize)];
    locData[(2*64)+localIdx] = ptrData[(i + localIdx) + (2*sourceSize)];
    locData[(3*64)+localIdx] = ptrData[(i + localIdx) + (3*sourceSize)];
    locData[(4*64)+localIdx] = ptrData[(i + localIdx) + (4*sourceSize)];

    locSource[localIdx] = ptrSource[i + localIdx];

    // Wait until all data is in shared memory
    __syncthreads();

    for (unsigned int k = 0; k < 64; k++) {
      curSupport = locSource[k];

      eval = locLevel_0 * locData[(0*64)+k];
      index_calc = eval - locIndex_0;
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel_1 * locData[(1*64)+k];
      index_calc = eval - locIndex_1;
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel_2 * locData[(2*64)+k];
      index_calc = eval - locIndex_2;
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel_3 * locData[(3*64)+k];
      index_calc = eval - locIndex_3;
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel_4 * locData[(4*64)+k];
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
                            unsigned int datasize,
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
    locData[d] = ptrData[(d*datasize)+globalIdx];
  }

  // Iterate over all grid points (fast ones, with cache)
  for (unsigned int j = 0; j < fastStorageSize; j += 64) {
#pragma unroll 5

    for (unsigned int d = 0; d < dims; d++) {
      locLevel[localIdx + (64 * d)] = ptrLevel[(j + localIdx) + (storageSize * d)];
      locIndex[localIdx + (64 * d)] = ptrIndex[(j + localIdx) + (storageSize * d)];
    }

    locAlpha[localIdx] = ptrAlpha[j + localIdx];

    // Wait until all needed data is in shared memory
    __syncthreads();

    for (unsigned int k = 0; k < 64; k++) {
      curSupport = locAlpha[k];

#pragma unroll 5

      for (unsigned int d = 0; d < dims; d++) {
        eval = locLevel[k+(64*d)] * locData[d];
        index_calc = eval - locIndex[k+(64*d)];
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
      eval = ptrLevel[m+(storageSize*d)] * locData[d];
      index_calc = eval - ptrIndex[m+(storageSize*d)];
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
                               unsigned int datasize,
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

  locData_0 = ptrData[(0*datasize)+globalIdx];
  locData_1 = ptrData[(1*datasize)+globalIdx];
  locData_2 = ptrData[(2*datasize)+globalIdx];
  locData_3 = ptrData[(3*datasize)+globalIdx];
  locData_4 = ptrData[(4*datasize)+globalIdx];

  // Iterate over all grid points (fast ones, with cache)
  for (unsigned int j = 0; j < fastStorageSize; j += 64) {
    locLevel[localIdx + (64*0)] = ptrLevel[(j + localIdx) + (storageSize*0)];
    locIndex[localIdx + (64*0)] = ptrIndex[(j + localIdx) + (storageSize*0)];
    locLevel[localIdx + (64*1)] = ptrLevel[(j + localIdx) + (storageSize*1)];
    locIndex[localIdx + (64*1)] = ptrIndex[(j + localIdx) + (storageSize*1)];
    locLevel[localIdx + (64*2)] = ptrLevel[(j + localIdx) + (storageSize*2)];
    locIndex[localIdx + (64*2)] = ptrIndex[(j + localIdx) + (storageSize*2)];
    locLevel[localIdx + (64*3)] = ptrLevel[(j + localIdx) + (storageSize*3)];
    locIndex[localIdx + (64*3)] = ptrIndex[(j + localIdx) + (storageSize*3)];
    locLevel[localIdx + (64*4)] = ptrLevel[(j + localIdx) + (storageSize*4)];
    locIndex[localIdx + (64*4)] = ptrIndex[(j + localIdx) + (storageSize*4)];

    locAlpha[localIdx] = ptrAlpha[j + localIdx];

    // Wait until all needed data is in shared memory
    __syncthreads();

    for (unsigned int k = 0; k < 64; k++) {
      curSupport = locAlpha[k];

      eval = locLevel[k + (0*64)] * locData_0;
      index_calc = eval - locIndex[k + (0*64)];
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel[k + (1*64)] * locData_1;
      index_calc = eval - locIndex[k + (1*64)];
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel[k + (2*64)] * locData_2;
      index_calc = eval - locIndex[k + (2*64)];
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel[k + (3*64)] * locData_3;
      index_calc = eval - locIndex[k + (3*64)];
      abs = fabsf(index_calc);
      last = 1.0f - abs;
      localSupport = fmaxf(last, 0.0f);
      curSupport *= localSupport;

      eval = locLevel[k + (4*64)] * locData_4;
      index_calc = eval - locIndex[k + (4*64)];
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

    eval = ptrLevel[m + (0*storageSize)] * locData_0;
    index_calc = eval - ptrIndex[m + (0*storageSize)];
    abs = fabsf(index_calc);
    last = 1.0f - abs;
    localSupport = fmaxf(last, 0.0f);
    curSupport *= localSupport;

    eval = ptrLevel[m + (1*storageSize)] * locData_1;
    index_calc = eval - ptrIndex[m + (1*storageSize)];
    abs = fabsf(index_calc);
    last = 1.0f - abs;
    localSupport = fmaxf(last, 0.0f);
    curSupport *= localSupport;

    eval = ptrLevel[m +(2*storageSize)] * locData_2;
    index_calc = eval - ptrIndex[m + (2*storageSize)];
    abs = fabsf(index_calc);
    last = 1.0f - abs;
    localSupport = fmaxf(last, 0.0f);
    curSupport *= localSupport;

    eval = ptrLevel[m + (3*storageSize)] * locData_3;
    index_calc = eval - ptrIndex[m + (3*storageSize)];
    abs = fabsf(index_calc);
    last = 1.0f - abs;
    localSupport = fmaxf(last, 0.0f);
    curSupport *= localSupport;

    eval = ptrLevel[m + (4*storageSize)] * locData_4;
    index_calc = eval - ptrIndex[m + (4*storageSize)];
    abs = fabsf(index_calc);
    last = 1.0f - abs;
    localSupport = fmaxf(last, 0.0f);
    curSupport *= localSupport;

    myResult += curSupport;
  }

  ptrResult[globalIdx] = myResult;
}

double multTransSPCUDA(float* ptrSource, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims) {
  // copy coefficients to GPU
  memcpy(host_datavec_sp, ptrSource, sourceSize * sizeof(float));
  cudaMemcpy(gpu_datavec_sp, ptrSource, sourceSize * sizeof(float), cudaMemcpyHostToDevice);
  int myStorageSize = ((int)storageSize) / 64;

  // Measure time
  sg::base::SGppStopwatch* mytimer = new sg::base::SGppStopwatch();
  mytimer->start();

  if (dims == 5) {
    multTransSP_CUDA_5d <<< myStorageSize, 64>>>(gpu_datavec_sp, gpu_dataset_sp, gpu_grid_level_sp, gpu_grid_index_sp, gpu_alpha_sp, (unsigned int)sourceSize, (unsigned int)gpu_full_storageSize,0);
  } else {
    multTransSP_CUDA <<< myStorageSize, 64>>>(gpu_datavec_sp, gpu_dataset_sp, gpu_grid_level_sp, gpu_grid_index_sp, gpu_alpha_sp, (unsigned int)sourceSize, (unsigned int)gpu_full_storageSize, (unsigned int)dims, 0);
  }

  // copy results back to host
  cudaMemcpy(host_alpha_sp, gpu_alpha_sp, storageSize * sizeof(float), cudaMemcpyDeviceToHost);
  memcpy(ptrGlobalResult, host_alpha_sp, storageSize * sizeof(float));

  double time = mytimer->stop();
  delete mytimer;

  return time;
}

double multSPCUDA(float* ptrAlpha, float* ptrResult, size_t result_size, size_t storageSize, size_t dims) {
  // copy coefficients to GPU
  memcpy(host_alpha_sp, ptrAlpha, storageSize * sizeof(float));
  cudaMemcpy(gpu_alpha_sp, host_alpha_sp, storageSize * sizeof(float), cudaMemcpyHostToDevice);

  unsigned int tmp = ((unsigned int)storageSize) / CUDA_BLOCK_SIZE_GPU;
  unsigned int fastStorageSize = tmp * CUDA_BLOCK_SIZE_GPU;

  // Measure time
  sg::base::SGppStopwatch* mytimer = new sg::base::SGppStopwatch();
  mytimer->start();

  if (dims == 5) {
    multSP_CUDA_5d <<< result_size / 64, 64 >>> (gpu_alpha_sp, gpu_dataset_sp, gpu_grid_level_sp, gpu_grid_index_sp, gpu_datavec_sp, fastStorageSize, (unsigned int)storageSize, (unsigned int)result_size, 0);
  } else {
    multSP_CUDA <<< result_size / 64, 64 >>> (gpu_alpha_sp, gpu_dataset_sp, gpu_grid_level_sp, gpu_grid_index_sp, gpu_datavec_sp, fastStorageSize, (unsigned int)storageSize, (unsigned int)dims, (unsigned int)result_size, 0);
  }

  // copy results back to host
  cudaMemcpy(host_datavec_sp, gpu_datavec_sp, result_size * sizeof(float), cudaMemcpyDeviceToHost);
  memcpy(ptrResult, host_datavec_sp, result_size * sizeof(float));

  double time = mytimer->stop();
  delete mytimer;

  return time;
}

void uploadGridSPCUDA(float* ptrLevel, float* ptrIndex, size_t storageSize, size_t dims) {
  size_t mem_size = storageSize * dims * sizeof(float);
  gpu_full_storageSize = storageSize;

  // allocate memory on GPU
  cudaMalloc((void**) &gpu_grid_level_sp, mem_size);
  cudaMalloc((void**) &gpu_grid_index_sp, mem_size);
  cudaMalloc((void**) &gpu_alpha_sp, storageSize * sizeof(float));
  cudaMallocHost((void**) &host_alpha_sp, storageSize * sizeof(float));
  cudaMallocHost((void**) &host_grid_level_sp, mem_size);
  cudaMallocHost((void**) &host_grid_index_sp, mem_size);

  memset(host_grid_level_sp, 0, mem_size);
  memset(host_grid_index_sp, 0, mem_size);

  // copy and transpose grid level and index to pinned memory
  //#pragma omp parallel for
  for (size_t i = 0; i < storageSize; i++) {
    for (size_t d = 0; d < dims; d++) {
      host_grid_level_sp[(d*storageSize)+i] = ptrLevel[(i*dims)+d];
      host_grid_index_sp[(d*storageSize)+i] = ptrIndex[(i*dims)+d];
    }
  } 

  // copy host memory to device
  cudaMemcpy(gpu_grid_level_sp, host_grid_level_sp, mem_size, cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_grid_index_sp, host_grid_index_sp, mem_size, cudaMemcpyHostToDevice);
}

void uploadDataSPCUDA(float* ptrData, size_t dataSize, size_t dims) {
  size_t mem_size = dataSize * dims * sizeof(float);

  // allocate memory on GPU
  cudaMalloc((void**) &gpu_dataset_sp, mem_size);
  cudaMalloc((void**) &gpu_datavec_sp, dataSize * sizeof(float));
  cudaMallocHost((void**) &host_datavec_sp, dataSize * sizeof(float));

  // copy host memory to device
  cudaMemcpy(gpu_dataset_sp, ptrData, mem_size, cudaMemcpyHostToDevice);
}

void deleteGridSPCUDA() {
  cudaFree(gpu_grid_level_sp);
  cudaFree(gpu_grid_index_sp);
  cudaFree(gpu_alpha_sp);
  cudaFreeHost(host_alpha_sp);
  cudaFreeHost(host_grid_level_sp);
  cudaFreeHost(host_grid_index_sp);
}

void deleteDataSPCUDA() {
  cudaFree(gpu_dataset_sp);
  cudaFree(gpu_datavec_sp);
  cudaFreeHost(host_datavec_sp);
}
