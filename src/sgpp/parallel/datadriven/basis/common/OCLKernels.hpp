/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OCLKERNELS_HPP
#define OCLKERNELS_HPP

#include <cstdlib>

#define MAX_OCL_DEVICE_COUNT 8
#define OCL_SGPP_LOCAL_WORKGROUP_SIZE 64

namespace sg
{
namespace parallel
{

class OCLKernels
{
private:
	bool isFirstTimeMultTransSP;
	bool isFirstTimeMultSP;
	bool isFirstTimeMultTransModSP;
	bool isFirstTimeMultModSP;

	bool isFirstTimeMultTransDP;
	bool isFirstTimeMultDP;
	bool isFirstTimeMultTransModDP;
	bool isFirstTimeMultModDP;

	bool isVeryFirstTimeDP;
	bool isVeryFirstTimeSP;

public:
	OCLKernels();

	~OCLKernels();

	void resetKernels();

	void resetData();

	double multTransOCL(double* ptrSource, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims, size_t gpu_partition);

	double multOCL(double* ptrAlpha, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrResult, size_t result_size, size_t storageSize, size_t dims, size_t gpu_partition);

	double multTransSPOCL(float* ptrSource, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims, size_t gpu_partition);

	double multSPOCL(float* ptrAlpha, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrResult, size_t result_size, size_t storageSize, size_t dims, size_t gpu_partition);

	double multTransModSPOCL(float* ptrSource, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims, size_t gpu_partition);

	double multModSPOCL(float* ptrAlpha, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrResult, size_t result_size, size_t storageSize, size_t dims, size_t gpu_partition);

	double multTransModOCL(double* ptrSource, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims, size_t gpu_partition);

	double multModOCL(double* ptrAlpha, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrResult, size_t result_size, size_t storageSize, size_t dims, size_t gpu_partition);
};

}

}

#endif /* OCLKERNELS_HPP */
