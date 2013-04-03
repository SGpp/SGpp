/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef CUDAKERNELS_HPP
#define CUDAKERNELS_HPP

#include <cstdlib>

#define CUDA_BLOCK_SIZE_GPU 64;
#define CUDA_PREFETCH_SIZE 64;

double multTransSPCUDA(float* ptrSource, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims);

double multSPCUDA(float* ptrAlpha, float* ptrResult, size_t result_size, size_t storageSize, size_t dims);

void uploadGridSPCUDA(float* ptrLevel, float* ptrIndex, size_t storageSize, size_t dims);

void uploadDataSPCUDA(float* ptrData, size_t dataSize, size_t dims);

void deleteGridSPCUDA();

void deleteDataSPCUDA();

#endif /* CUDAKERNELS_HPP */
