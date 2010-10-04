/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef ARBBKERNELS_HPP
#define ARBBKERNELS_HPP

#include <cstdlib>

namespace sg
{

class ArBBKernels
{
private:
	bool isMultTransSPfirst;
	bool isMultSPfirst;
	bool isMultTransfirst;
	bool isMultfirst;

public:
	ArBBKernels();

	~ArBBKernels();

	double multArBB(double* ptrSource, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims);

	double multTransArBB(double* ptrAlpha, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrResult, size_t result_size, size_t storageSize, size_t dims);

	double multSPArBB(float* ptrSource, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims);

	double multTransSPArBB(float* ptrAlpha, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrResult, size_t result_size, size_t storageSize, size_t dims);

	void resetKernels();
};

}

#endif /* ARBBKERNELS_HPP */
