/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/noboundary/operation/datadriven/ArBBKernels.hpp"

#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <iostream>
#include <string>

#include <arbb.hpp>

//#define OLD_ARBB_IMPL

namespace sg
{

arbb::dense<arbb::f32, 2> ArBB_DataSP;
arbb::dense<arbb::f32, 2> ArBB_LevelSP;
arbb::dense<arbb::f32, 2> ArBB_IndexSP;

arbb::dense<arbb::f64, 2> ArBB_Data;
arbb::dense<arbb::f64, 2> ArBB_Level;
arbb::dense<arbb::f64, 2> ArBB_Index;

#ifdef OLD_ARBB_IMPL
template <typename fp_Type>
void arbb_evalGridPoint(const arbb::dense<fp_Type>& dataPoint, const arbb::dense<fp_Type>& levelPoint, const arbb::dense<fp_Type>& indexPoint, fp_Type& result)
{
	fp_Type eval;
	fp_Type index_calc;
	fp_Type abs;
	fp_Type last;
	fp_Type localSupport;
	fp_Type one = 1.0;
	fp_Type zero = 0.0;
	arbb::usize dims = dataPoint.num_cols();

	_for (arbb::usize i = 0, i < dims, i++)
	{
		eval = levelPoint[i]*dataPoint[i];
		index_calc = eval - indexPoint[i];
		abs = arbb::abs(index_calc);
		last = one - abs;
		localSupport = arbb::max(last, zero);
		result *= localSupport;
	} _end_for;
}
#endif

template <typename fp_Type>
void arbb_evalGridPoint_oneDim(const arbb::dense<fp_Type>& dataPointsDim, const fp_Type& levelPoint, const fp_Type& indexPoint, arbb::dense<fp_Type>& result)
{
	arbb::dense<fp_Type> zero = arbb::fill(static_cast<fp_Type>(0), result.length());
	arbb::dense<fp_Type> one = arbb::fill(static_cast<fp_Type>(1), result.length());
	arbb::dense<fp_Type> index = arbb::fill(indexPoint, result.length());

	result = levelPoint * dataPointsDim;
	result = result - index;
	result = arbb::abs(result);
	result = one - result;
	result = arbb::max(result, zero);
}

template <typename fp_Type>
void arbb_evalTransGridPoint_oneDim(const fp_Type& dataPointDim, const arbb::dense<fp_Type>& level, const arbb::dense<fp_Type>& index, arbb::dense<fp_Type>& result)
{
	arbb::dense<fp_Type> zero = arbb::fill(static_cast<fp_Type>(0), result.length());
	arbb::dense<fp_Type> one = arbb::fill(static_cast<fp_Type>(1), result.length());

	result = level * dataPointDim;
	result = result - index;
	result = arbb::abs(result);
	result = one - result;
	result = arbb::max(result, zero);
}

template <typename fp_Type>
void arbb_mult(const arbb::dense<fp_Type, 2>& Data, const arbb::dense<fp_Type, 2>& Level, const arbb::dense<fp_Type, 2>& Index, const arbb::dense<fp_Type>& source, arbb::dense<fp_Type>& result)
{
#ifdef OLD_ARBB_IMPL
	arbb::usize source_size = Data.num_rows();
	arbb::usize storage_size = Level.num_rows();

	_for (arbb::usize i = 0, i < source_size, i++)
	{
		_for (arbb::usize j = 0, j < storage_size, j++)
		{
			fp_Type curResult = source[i];
			fp_Type lastResult = result[j];

			arbb_evalGridPoint(Data.row(i), Level.row(j), Index.row(j), curResult);

			lastResult += curResult;
			result[j] = lastResult;
		} _end_for;
	} _end_for;
#else
	arbb::usize dims = Level.num_cols();
	arbb::usize source_size = Data.num_rows();
	arbb::usize storage_size = Level.num_rows();

	_for (arbb::usize i = 0, i < source_size, i++)
	{
		arbb::dense<fp_Type> curDataPoint = Data.row(i);
		fp_Type a = source[i];

		arbb::dense<fp_Type> temp_result = arbb::fill(a, storage_size);
		arbb::dense<fp_Type> temp = arbb::fill(static_cast<fp_Type>(0), storage_size);

		_for (arbb::usize d = 0, d < dims, d++)
		{
			fp_Type dataDim = curDataPoint[d];

			arbb_evalTransGridPoint_oneDim(dataDim, Level.col(d), Index.col(d), temp);

			temp_result *= temp;
		} _end_for;

		result += temp_result;
	} _end_for;

#endif
}

template <typename fp_Type>
void arbb_multTrans(const arbb::dense<fp_Type, 2>& Data, const arbb::dense<fp_Type, 2>& Level, const arbb::dense<fp_Type, 2>& Index, const arbb::dense<fp_Type>& alpha, arbb::dense<fp_Type>& result)
{
#ifdef OLD_ARBB_IMPL
	arbb::usize result_size = Data.num_rows();
	arbb::usize storage_size = Level.num_rows();

	_for (arbb::usize i = 0, i < result_size, i++)
	{
		fp_Type lastResult = result[i];

		_for (arbb::usize j = 0, j < storage_size, j++)
		{
			fp_Type curResult = alpha[j];
			arbb_evalGridPoint(Data.row(i), Level.row(j), Index.row(j), curResult);
			lastResult += curResult;
		} _end_for;

		result[i] = lastResult;
	} _end_for;
#else
	arbb::usize dims = Level.num_cols();
	arbb::usize result_size = Data.num_rows();
	arbb::usize storage_size = Level.num_rows();

	_for (arbb::usize j = 0, j < storage_size, j++)
	{
		arbb::dense<fp_Type> curLevel = Level.row(j);
		arbb::dense<fp_Type> curIndex = Index.row(j);
		fp_Type a = alpha[j];

		arbb::dense<fp_Type> temp_result = arbb::fill(a, result_size);
		arbb::dense<fp_Type> temp = arbb::fill(static_cast<fp_Type>(0), result_size);

		_for (arbb::usize d = 0, d < dims, d++)
		{
			fp_Type l = curLevel[d];
			fp_Type i = curIndex[d];

			arbb_evalGridPoint_oneDim(Data.col(d), l, i, temp);

			temp_result *= temp;
		} _end_for;

		result += temp_result;
	} _end_for;
#endif
}

ArBBKernels::ArBBKernels()
{
	isMultTransSPfirst = true;
	isMultSPfirst = true;

	isMultTransfirst = true;
	isMultfirst = true;
}

ArBBKernels::~ArBBKernels()
{
}

double ArBBKernels::multArBB(double* ptrSource, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims)
{
	double time = 0.0;

	try
	{
		arbb::dense<arbb::f64> ArBB_result;
		arbb::dense<arbb::f64> ArBB_source;

		if (isMultTransfirst && isMultfirst)
		{
			arbb::bind(ArBB_Data, ptrData, dims, sourceSize);
			arbb::bind(ArBB_Level, ptrLevel, dims, storageSize);
			arbb::bind(ArBB_Index, ptrIndex, dims, storageSize);
			isMultfirst = false;
		}

		arbb::bind(ArBB_result, ptrGlobalResult, storageSize);
		arbb::bind(ArBB_source, ptrSource, sourceSize);

		arbb::call(&(arbb_mult<arbb::f64>))(ArBB_Data, ArBB_Level, ArBB_Index, ArBB_source, ArBB_result);
	}
	catch (const std::exception& e)
	{
		std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
	}

//   	for (size_t i = 0; i < sourceSize; i++)
//	{
//		for (size_t j = 0; j < storageSize; j++)
//		{
//			double curSupport = ptrSource[i];
//
//			for (size_t d = 0; d < dims; d++)
//			{
//				double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d]));
//				double index_calc = eval - (ptrIndex[(j*dims)+d]);
//				double abs = fabs(index_calc);
//				double last = 1.0 - abs;
//				double localSupport = std::max<double>(last, 0.0);
//				curSupport *= localSupport;
//			}
//
//			ptrGlobalResult[j] += curSupport;
//		}
//	}

	return time;
}

double ArBBKernels::multTransArBB(double* ptrAlpha, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrResult, size_t result_size, size_t storageSize, size_t dims)
{
	double time = 0.0;

	try
	{
		arbb::dense<arbb::f64> ArBB_result;
		arbb::dense<arbb::f64> ArBB_alpha;

		if (isMultTransfirst && isMultfirst)
		{
			arbb::bind(ArBB_Data, ptrData, dims, result_size);
			arbb::bind(ArBB_Level, ptrLevel, dims, storageSize);
			arbb::bind(ArBB_Index, ptrIndex, dims, storageSize);
			isMultTransfirst = false;
		}

		arbb::bind(ArBB_result, ptrResult, result_size);
		arbb::bind(ArBB_alpha, ptrAlpha, storageSize);

		arbb::call(&(arbb_multTrans<arbb::f64>))(ArBB_Data, ArBB_Level, ArBB_Index, ArBB_alpha, ArBB_result);
	}
	catch (const std::exception& e)
	{
		std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
	}

//	for (size_t i = 0; i < result_size; i++)
//	{
//		for (size_t j = 0; j < storageSize; j++)
//		{
//			double curSupport = ptrAlpha[j];
//
//			for (size_t d = 0; d < dims; d++)
//			{
//				double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d]));
//				double index_calc = eval - (ptrIndex[(j*dims)+d]);
//				double abs = fabs(index_calc);
//				double last = 1.0 - abs;
//				double localSupport = std::max<double>(last, 0.0);
//				curSupport *= localSupport;
//			}
//
//			ptrResult[i] += curSupport;
//		}
//	}

	return time;
}

double ArBBKernels::multSPArBB(float* ptrSource, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims)
{
	double time = 0.0;

	try
	{
		arbb::dense<arbb::f32> ArBB_result;
		arbb::dense<arbb::f32> ArBB_source;

		if (isMultTransSPfirst && isMultSPfirst)
		{
			arbb::bind(ArBB_DataSP, ptrData, dims, sourceSize);
			arbb::bind(ArBB_LevelSP, ptrLevel, dims, storageSize);
			arbb::bind(ArBB_IndexSP, ptrIndex, dims, storageSize);
			isMultSPfirst = false;
		}

		arbb::bind(ArBB_result, ptrGlobalResult, storageSize);
		arbb::bind(ArBB_source, ptrSource, sourceSize);

		arbb::call(&(arbb_mult<arbb::f32>))(ArBB_DataSP, ArBB_LevelSP, ArBB_IndexSP, ArBB_source, ArBB_result);
	}
	catch (const std::exception& e)
	{
		std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
	}

//   	for (size_t i = 0; i < sourceSize; i++)
//	{
//		for (size_t j = 0; j < storageSize; j++)
//		{
//			float curSupport = ptrSource[i];
//
//			for (size_t d = 0; d < dims; d++)
//			{
//				float eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d]));
//				float index_calc = eval - (ptrIndex[(j*dims)+d]);
//				float abs = fabs(index_calc);
//				float last = 1.0f - abs;
//				float localSupport = std::max<float>(last, 0.0f);
//				curSupport *= localSupport;
//			}
//
//			ptrGlobalResult[j] += curSupport;
//		}
//	}

	return time;
}

double ArBBKernels::multTransSPArBB(float* ptrAlpha, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrResult, size_t result_size, size_t storageSize, size_t dims)
{
	double time = 0.0;

	try
	{
		arbb::dense<arbb::f32> ArBB_result;
		arbb::dense<arbb::f32> ArBB_alpha;

		if (isMultTransSPfirst && isMultSPfirst)
		{
			arbb::bind(ArBB_DataSP, ptrData, dims, result_size);
			arbb::bind(ArBB_LevelSP, ptrLevel, dims, storageSize);
			arbb::bind(ArBB_IndexSP, ptrIndex, dims, storageSize);
			isMultTransSPfirst = false;
		}

		arbb::bind(ArBB_result, ptrResult, result_size);
		arbb::bind(ArBB_alpha, ptrAlpha, storageSize);

		arbb::call(&(arbb_multTrans<arbb::f32>))(ArBB_DataSP, ArBB_LevelSP, ArBB_IndexSP, ArBB_alpha, ArBB_result);
	}
	catch (const std::exception& e)
	{
		std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
	}

//	for (size_t i = 0; i < result_size; i++)
//	{
//		for (size_t j = 0; j < storageSize; j++)
//		{
//			float curSupport = ptrAlpha[j];
//
//			for (size_t d = 0; d < dims; d++)
//			{
//				float eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d]));
//				float index_calc = eval - (ptrIndex[(j*dims)+d]);
//				float abs = fabs(index_calc);
//				float last = 1.0f - abs;
//				float localSupport = std::max<double>(last, 0.0f);
//				curSupport *= localSupport;
//			}
//
//			ptrResult[i] += curSupport;
//		}
//	}

	return time;
}

void ArBBKernels::resetKernels()
{
	isMultTransSPfirst = true;
	isMultSPfirst = true;

	isMultTransfirst = true;
	isMultfirst = true;
}

}
