/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/noboundary/operation/datadriven/OCLKernels.hpp"

#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <iostream>

namespace sg
{

double multOCL(double* ptrSource, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims)
{
	double time = 0.0;

   	for (size_t i = 0; i < sourceSize; i++)
	{
		for (size_t j = 0; j < storageSize; j++)
		{
			double curSupport = ptrSource[i];

			for (size_t d = 0; d < dims; d++)
			{
				double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*sourceSize)+i]));
				double index_calc = eval - (ptrIndex[(j*dims)+d]);
				double abs = fabs(index_calc);
				double last = 1.0 - abs;
				double localSupport = std::max<double>(last, 0.0);
				curSupport *= localSupport;
			}

			ptrGlobalResult[j] += curSupport;
		}
	}

	return time;
}

double multTransOCL(double* ptrAlpha, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrResult, size_t result_size, size_t storageSize, size_t dims)
{
	double time = 0.0;

	for (size_t i = 0; i < result_size; i++)
	{
		for (size_t j = 0; j < storageSize; j++)
		{
			double curSupport = ptrAlpha[j];

			for (size_t d = 0; d < dims; d++)
			{
				double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*result_size)+i]));
				double index_calc = eval - (ptrIndex[(j*dims)+d]);
				double abs = fabs(index_calc);
				double last = 1.0 - abs;
				double localSupport = std::max<double>(last, 0.0);
				curSupport *= localSupport;
			}

			ptrResult[i] += curSupport;
		}
	}

	return time;
}

double multSPOCL(float* ptrSource, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims)
{
	double time = 0.0;

   	for (size_t i = 0; i < sourceSize; i++)
	{
		for (size_t j = 0; j < storageSize; j++)
		{
			float curSupport = ptrSource[i];

			for (size_t d = 0; d < dims; d++)
			{
				float eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*sourceSize)+i]));
				float index_calc = eval - (ptrIndex[(j*dims)+d]);
				float abs = fabs(index_calc);
				float last = 1.0 - abs;
				float localSupport = std::max<float>(last, 0.0);
				curSupport *= localSupport;
			}

			ptrGlobalResult[j] += curSupport;
		}
	}

	return time;
}

double multTransSPOCL(float* ptrAlpha, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrResult, size_t result_size, size_t storageSize, size_t dims)
{
	double time = 0.0;

	for (size_t i = 0; i < result_size; i++)
	{
		for (size_t j = 0; j < storageSize; j++)
		{
			float curSupport = ptrAlpha[j];

			for (size_t d = 0; d < dims; d++)
			{
				float eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*result_size)+i]));
				float index_calc = eval - (ptrIndex[(j*dims)+d]);
				float abs = fabs(index_calc);
				float last = 1.0 - abs;
				float localSupport = std::max<float>(last, 0.0);
				curSupport *= localSupport;
			}

			ptrResult[i] += curSupport;
		}
	}

	return time;
}

}
