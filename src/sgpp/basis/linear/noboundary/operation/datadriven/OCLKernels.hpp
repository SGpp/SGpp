/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OCLKERNELS_HPP
#define OCLKERNELS_HPP

#include <cstdlib>

#define OCL_MULT_N_DATAPREFETCH_BLOCKSIZE 32 // must divide vec-length with no reminder
#define OCL_DATAPREFETCH_SIZE 64 // must divide vec-length with no reminder

namespace sg
{

class OCLKernels
{
private:
	bool isFirstTimeMultTransSP;
	bool isFirstTimeMultSP;
	bool isVeryFirstTimeSP;

public:
	OCLKernels();

	~OCLKernels();

	void resetKernels();

	void resetData();

	double multOCL(double* ptrSource, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims);

	double multTransOCL(double* ptrAlpha, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrResult, size_t result_size, size_t storageSize, size_t dims);

	double multSPOCL(float* ptrSource, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims);

	double multTransSPOCL(float* ptrAlpha, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrResult, size_t result_size, size_t storageSize, size_t dims);
};

}

#endif /* OCLKERNELS_HPP */
