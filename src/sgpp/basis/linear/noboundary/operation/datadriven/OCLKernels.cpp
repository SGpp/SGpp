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
#include <string>
#include <sstream>

#define OCL_MULT_BLOCKSIZE 32
#define OCL_MULTTRANS_BLOCKSIZE 64

// include OpenCL
#include "CL/cl.h"

namespace sg
{

cl_int err;
cl_platform_id platform_id;
cl_device_id device_id;
cl_uint num_platforms;
cl_uint num_devices;
cl_context context;
cl_command_queue command_queue;

cl_mem clDataSP, clLevelSP, clIndexSP;

cl_kernel kernel_multTransSP;
cl_program program_multTransSP;

cl_kernel kernel_multSP;
cl_program program_multSP;

OCLKernels::OCLKernels()
{
	// Check for OCL Platforms
	err = clGetPlatformIDs(1, &platform_id, &num_platforms);
	if (err != CL_SUCCESS)
	{
		std::cout << "Unable to get Platform ID. Error Code: " << err << std::endl;
	}

	// Find out how many devices there are
	err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &num_devices);
    if (err != CL_SUCCESS)
    {
    	std::cout << "Unable to get Device ID. Error Code: " << err << std::endl;
    }
    else if (num_devices == 0)
    {
    	std::cout << "NO GPU OpenCL devices have been found!" << std::endl;
    }

    // Create GPU context
    context = clCreateContext(0, num_devices, &device_id, NULL, NULL, &err);
    if (err != CL_SUCCESS)
	{
    	std::cout << "Failed to create OpenCL context! Error Code: " << err << std::endl;
	}

    // Creating the command queue
    command_queue = clCreateCommandQueue(context, device_id, 0, &err);
    if (err != CL_SUCCESS)
	{
    	std::cout << "Failed to create command queue! Error Code: " << err << std::endl;
	}

    isFirstTimeMultTransSP = true;
    isFirstTimeMultSP = true;
}

OCLKernels::~OCLKernels()
{
	if (!isFirstTimeMultTransSP)
	{
	    clReleaseProgram(program_multTransSP);
	    clReleaseKernel(kernel_multTransSP);
	}

	if (!isFirstTimeMultSP)
	{
		clReleaseMemObject(clDataSP);
		clReleaseMemObject(clLevelSP);
		clReleaseMemObject(clIndexSP);
	    clReleaseProgram(program_multSP);
	    clReleaseKernel(kernel_multSP);
	}

    clReleaseCommandQueue(command_queue);
    clReleaseContext(context);
}

double OCLKernels::multOCL(double* ptrSource, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims)
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

double OCLKernels::multTransOCL(double* ptrAlpha, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrResult, size_t result_size, size_t storageSize, size_t dims)
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

double OCLKernels::multSPOCL(float* ptrSource, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims)
{
	double time = 0.0;

	if (isFirstTimeMultSP)
	{
		std::stringstream stream_program_src;

		stream_program_src << "__kernel void multSPOCL(__global float* ptrSource," << std::endl;
		stream_program_src << "						__global float* ptrData," << std::endl;
		stream_program_src << "						__global float* ptrLevel," << std::endl;
		stream_program_src << "						__global float* ptrIndex," << std::endl;
		stream_program_src << "						__global float* ptrResult," << std::endl;
		stream_program_src << "						uint sourceSize)" << std::endl;
		stream_program_src << "{" << std::endl;
		stream_program_src << "	int globalIdx = get_global_id(0);" << std::endl << std::endl;
		stream_program_src << "	float eval, index_calc, abs, last, localSupport;" << std::endl << std::endl;
		stream_program_src << "	float myResult = 0.0f;" << std::endl << std::endl;
		stream_program_src << "	// Iterate over all grid points" << std::endl;
		stream_program_src << "	for(int i = 0; i < sourceSize; i++)" << std::endl;
		stream_program_src << "	{" << std::endl;
		stream_program_src << "		float curSupport = ptrSource[i];" << std::endl << std::endl;
		for (size_t d = 0; d < dims; d++)
		{
			//stream_program_src << "			prefetch(&ptrData[(" << d << "*" << sourceSize << ")+(i+20)], 1);" << std::endl;
			stream_program_src << "			eval = ((ptrLevel[(globalIdx*" << dims << ")+" << d << "]) * (ptrData[(i*" << dims << ")+" << d << "]));" << std::endl;
			stream_program_src << "			index_calc = eval - (ptrIndex[(globalIdx*" << dims << ")+" << d << "]);" << std::endl;
			stream_program_src << "			abs = fabs(index_calc);" << std::endl;
			stream_program_src << "			last = 1.0f - abs;" << std::endl;
			stream_program_src << "			localSupport = fmax(last, 0.0f);" << std::endl;
			stream_program_src << "			curSupport *= localSupport;" << std::endl;
		}
		stream_program_src << std::endl << "		myResult += curSupport;" << std::endl;
		stream_program_src << "	}" << std::endl;
		stream_program_src << "	ptrResult[globalIdx] = myResult;" << std::endl;
		stream_program_src << "}" << std::endl;

		std::string program_src = stream_program_src.str();

		//std::cout << program_src << std::endl;

	    // setting the program
	    const char* kernel_src = program_src.c_str();
	    program_multSP = clCreateProgramWithSource(context, 1, &kernel_src, NULL, &err);
	    if (err != CL_SUCCESS)
		{
	    	std::cout << "Failed to create program! Error Code: " << err << std::endl;
	    	return 0.0;
		}

	    // compiling the program
	    err = clBuildProgram(program_multSP, 0, NULL, "-cl-finite-math-only -cl-strict-aliasing -cl-fast-relaxed-math -cl-single-precision-constant", NULL, NULL);
	    //err = clBuildProgram(program_multSP, 0, NULL, "-cl-opt-disable", NULL, NULL);
	    //err = clBuildProgram(program_multSP, 0, NULL, NULL, NULL, NULL);
	    if (err != CL_SUCCESS)
	    {
	    	std::cout << "OpenCL Build Error. Error Code: " << err << std::endl;

	    	size_t len;
	    	char buffer[2048];

	    	// get the build log
	    	clGetProgramBuildInfo(program_multSP, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);

	    	std::cout << "--- Build Log ---" << std::endl << buffer << std::endl;
	    	return 0.0;
	    }

	    // creating the kernel
	    kernel_multSP = clCreateKernel(program_multSP, "multSPOCL", &err);
	    if (err != CL_SUCCESS)
		{
	    	std::cout << "Failed to create kernel! Error Code: " << err << std::endl;
	    	return 0.0;
		}
	}

	if (isFirstTimeMultSP && isFirstTimeMultTransSP)
	{
		clDataSP = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*dims*sourceSize, ptrData, NULL);
		clLevelSP = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*dims*storageSize, ptrLevel, NULL);
		clIndexSP = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*dims*storageSize, ptrIndex, NULL);
	}

	cl_mem clSource, clResult;
	clSource = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*sourceSize, ptrSource, NULL);
    clResult = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float)*storageSize, NULL, NULL);

    cl_uint clSourceSize = (cl_uint)sourceSize;
	// set kernel arguments
	if ( clSetKernelArg(kernel_multSP, 0, sizeof(cl_mem), &clSource) ||
			clSetKernelArg(kernel_multSP, 1, sizeof(cl_mem), &clDataSP) ||
			clSetKernelArg(kernel_multSP, 2, sizeof(cl_mem), &clLevelSP) ||
			clSetKernelArg(kernel_multSP, 3, sizeof(cl_mem), &clIndexSP) ||
			clSetKernelArg(kernel_multSP, 4, sizeof(cl_mem), &clResult) ||
			clSetKernelArg(kernel_multSP, 5, sizeof(cl_uint), &clSourceSize) != CL_SUCCESS)
	{
		std::cout << "Failed to create kernel Args!" << std::endl;
		return 0.0;
	}

    // determine best fit devideable by OCL_MULT_BLOCKSIZE
	size_t numWGs = storageSize/OCL_MULT_BLOCKSIZE;
    size_t global = numWGs*OCL_MULT_BLOCKSIZE;
    if (global == 0)
    {
    	global = storageSize;
    }
    //size_t global = storageSize;
    // enqueue kernel
    err = clEnqueueNDRangeKernel(command_queue, kernel_multSP, 1, NULL, &global, NULL, 0, NULL, NULL);
    if (err != CL_SUCCESS)
	{
    	std::cout << "Failed to enqueue kernel command! Error Code: " << err << std::endl;
    	return 0.0;
	}

    // wait for command to finish
    clFinish(command_queue);

    // read data back
    err = clEnqueueReadBuffer(command_queue, clResult, CL_TRUE, 0, sizeof(float)*global, ptrGlobalResult, 0, NULL, NULL);
    if (err != CL_SUCCESS)
	{
    	std::cout << "Failed to enqueue read buffer command (mult)! Error Code: " << err << std::endl;
    	return 0.0;
	}

    // clean up
    clReleaseMemObject(clSource);
    clReleaseMemObject(clResult);

    isFirstTimeMultSP = false;

    // do the rest...
	for (size_t j = global; j < storageSize; j++)
	{
		ptrGlobalResult[j] = 0.0f;

		for (size_t i = 0; i < sourceSize; i++)
		{
			float curSupport = ptrSource[i];

			for (size_t d = 0; d < dims; d++)
			{
				float eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d]));
				float index_calc = eval - (ptrIndex[(j*dims)+d]);
				float abs = fabs(index_calc);
				float last = 1.0f - abs;
				float localSupport = std::max<float>(last, 0.0f);
				curSupport *= localSupport;
			}

			ptrGlobalResult[j] += curSupport;
		}
	}

	return time;
}

double OCLKernels::multTransSPOCL(float* ptrAlpha, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrResult, size_t result_size, size_t storageSize, size_t dims)
{
	double time = 0.0;

	if (isFirstTimeMultTransSP)
	{
		std::stringstream stream_program_src;

		stream_program_src << "__kernel void multTransSPOCL(__global float* ptrAlpha," << std::endl;
		stream_program_src << "						__global float* ptrData," << std::endl;
		stream_program_src << "						__global float* ptrLevel," << std::endl;
		stream_program_src << "						__global float* ptrIndex," << std::endl;
		stream_program_src << "						__global float* ptrResult," << std::endl;
		stream_program_src << "						uint fastStorageSize," << std::endl;
		stream_program_src << "						uint storageSize)" << std::endl;
		stream_program_src << "{" << std::endl;
		stream_program_src << "	int globalIdx = get_global_id(0);" << std::endl;
		stream_program_src << "	int localIdx = get_local_id(0);" << std::endl;
		stream_program_src << std::endl;
		stream_program_src << "	__local float locLevel[" << dims * OCL_MULTTRANS_BLOCKSIZE << "];" << std::endl;
		stream_program_src << "	__local float locIndex[" << dims * OCL_MULTTRANS_BLOCKSIZE << "];" << std::endl;
		stream_program_src << std::endl;
		stream_program_src << "	float eval, index_calc, abs, last, localSupport, curSupport;" << std::endl << std::endl;
		stream_program_src << "	float myResult = 0.0f;" << std::endl << std::endl;
		stream_program_src << "	// Create registers for the data" << std::endl;
		for (size_t d = 0; d < dims; d++)
		{
			stream_program_src << "	float data_" << d << " = ptrData[(globalIdx*" << dims << ")+" << d << "];" << std::endl;
		}
		stream_program_src << std::endl;
		stream_program_src << "	// Iterate over all grid points (fast ones, with cache)" << std::endl;
		stream_program_src << "	for(int j = 0; j < fastStorageSize; j+=" << OCL_MULTTRANS_BLOCKSIZE << ")" << std::endl;
		stream_program_src << "	{" << std::endl;
		for (size_t d = 0; d < dims; d++)
		{
			stream_program_src << "		locLevel[(localIdx*" << dims << ")+"<< d << "] = ptrLevel[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
			stream_program_src << "		locIndex[(localIdx*" << dims << ")+"<< d << "] = ptrIndex[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
		}
		stream_program_src << "		barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
		stream_program_src << std::endl;
		stream_program_src << "		for(int k = 0; k < " << OCL_MULTTRANS_BLOCKSIZE << "; k++)" << std::endl;
		stream_program_src << "		{" << std::endl;
		stream_program_src << "			curSupport = ptrAlpha[j+k];" << std::endl << std::endl;
		for (size_t d = 0; d < dims; d++)
		{
			stream_program_src << "			eval = ((locLevel[(k*" << dims << ")+" << d << "]) * (data_" << d << "));" << std::endl;
			stream_program_src << "			index_calc = eval - (locIndex[(k*" << dims << ")+" << d << "]);" << std::endl;
			stream_program_src << "			abs = fabs(index_calc);" << std::endl;
			stream_program_src << "			last = 1.0f - abs;" << std::endl;
			stream_program_src << "			localSupport = fmax(last, 0.0f);" << std::endl;
			stream_program_src << "			curSupport *= localSupport;" << std::endl << std::endl;
		}
		stream_program_src << "			myResult += curSupport;" << std::endl;
		stream_program_src << "		}" << std::endl;
		stream_program_src << std::endl;
		stream_program_src << "		barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
		stream_program_src << "	}" << std::endl;
		stream_program_src << std::endl;
		stream_program_src << "	// Iterate over all grid points (slow ones, without cache)" << std::endl;
		stream_program_src << "	for(int m = fastStorageSize; m < storageSize; m++)" << std::endl;
		stream_program_src << "	{" << std::endl;
		stream_program_src << "		curSupport = ptrAlpha[m];" << std::endl << std::endl;
		for (size_t d = 0; d < dims; d++)
		{
			stream_program_src << "		eval = ((ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << "));" << std::endl;
			stream_program_src << "		index_calc = eval - (ptrIndex[(m*" << dims << ")+" << d << "]);" << std::endl;
			stream_program_src << "		abs = fabs(index_calc);" << std::endl;
			stream_program_src << "		last = 1.0f - abs;" << std::endl;
			stream_program_src << "		localSupport = fmax(last, 0.0f);" << std::endl;
			stream_program_src << "		curSupport *= localSupport;" << std::endl << std::endl;
		}
		stream_program_src << "		myResult += curSupport;" << std::endl;
		stream_program_src << "	}" << std::endl;
		stream_program_src << std::endl;
		stream_program_src << "	ptrResult[globalIdx] = myResult;" << std::endl;
		stream_program_src << "}" << std::endl;

		std::string program_src = stream_program_src.str();

		//std::cout << program_src << std::endl;

	    // setting the program
	    const char* kernel_src = program_src.c_str();
	    program_multTransSP = clCreateProgramWithSource(context, 1, &kernel_src, NULL, &err);
	    if (err != CL_SUCCESS)
		{
	    	std::cout << "Failed to create program! Error Code: " << err << std::endl;
	    	return 0.0;
		}

	    // compiling the program
	    err = clBuildProgram(program_multTransSP, 0, NULL,  "-cl-finite-math-only -cl-strict-aliasing -cl-fast-relaxed-math -cl-single-precision-constant", NULL, NULL);
	    //err = clBuildProgram(program_multTransSP, 0, NULL, "-cl-opt-disable", NULL, NULL);
	    //err = clBuildProgram(program_multTransSP, 0, NULL, NULL, NULL, NULL);
	    if (err != CL_SUCCESS)
	    {
	    	std::cout << "OpenCL Build Error. Error Code: " << err << std::endl;

	    	size_t len;
	    	char buffer[2048];

	    	// get the build log
	    	clGetProgramBuildInfo(program_multTransSP, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);

	    	std::cout << "--- Build Log ---" << std::endl << buffer << std::endl;
	    	return 0.0;
	    }

	    // creating the kernel
	    kernel_multTransSP = clCreateKernel(program_multTransSP, "multTransSPOCL", &err);
	    if (err != CL_SUCCESS)
		{
	    	std::cout << "Failed to create kernel! Error Code: " << err << std::endl;
	    	return 0.0;
		}
	}

	if (isFirstTimeMultSP && isFirstTimeMultTransSP)
	{
		clDataSP = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*dims*result_size, ptrData, NULL);
		clLevelSP = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*dims*storageSize, ptrLevel, NULL);
		clIndexSP = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*dims*storageSize, ptrIndex, NULL);
	}

	cl_mem clAlpha, clResult;
    clAlpha = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*storageSize, ptrAlpha, NULL);
    clResult = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float)*result_size, NULL, NULL);

    size_t oclStorageSize = storageSize/OCL_MULTTRANS_BLOCKSIZE;
    oclStorageSize *= OCL_MULTTRANS_BLOCKSIZE;
    cl_uint clFastStorageSize = (cl_uint)(oclStorageSize);
    cl_uint clStorageSize = (cl_uint)(storageSize);

	// set kernel arguments
	if ( clSetKernelArg(kernel_multTransSP, 0, sizeof(cl_mem), &clAlpha) ||
			clSetKernelArg(kernel_multTransSP, 1, sizeof(cl_mem), &clDataSP) ||
			clSetKernelArg(kernel_multTransSP, 2, sizeof(cl_mem), &clLevelSP) ||
			clSetKernelArg(kernel_multTransSP, 3, sizeof(cl_mem), &clIndexSP) ||
			clSetKernelArg(kernel_multTransSP, 4, sizeof(cl_mem), &clResult) ||
			clSetKernelArg(kernel_multTransSP, 5, sizeof(cl_uint), &clFastStorageSize) ||
			clSetKernelArg(kernel_multTransSP, 6, sizeof(cl_uint), &clStorageSize) != CL_SUCCESS)
	{
		std::cout << "Failed to create kernel Args!" << std::endl;
		return 0.0;
	}

    size_t global = result_size;
    size_t local = OCL_MULTTRANS_BLOCKSIZE;
    // enqueue kernel
    err = clEnqueueNDRangeKernel(command_queue, kernel_multTransSP, 1, NULL, &global, &local, 0, NULL, NULL);
    if (err != CL_SUCCESS)
	{
    	std::cout << "Failed to enqueue kernel command! Error Code: " << err << std::endl;
    	return 0.0;
	}

    // wait for command to finish
    clFinish(command_queue);

    // read data back
    err = clEnqueueReadBuffer(command_queue, clResult, CL_TRUE, 0, sizeof(float)*result_size, ptrResult, 0, NULL, NULL);
    if (err != CL_SUCCESS)
	{
    	std::cout << "Failed to enqueue read buffer command (multTrans)! Error Code: " << err << std::endl;
    	return 0.0;
	}

    // clean up
    clReleaseMemObject(clAlpha);
    clReleaseMemObject(clResult);

    isFirstTimeMultTransSP = false;

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
//				float localSupport = std::max<float>(last, 0.0f);
//				curSupport *= localSupport;
//			}
//
//			ptrResult[i] += curSupport;
//		}
//	}

	return time;
}

void OCLKernels::resetKernels()
{
	if (!isFirstTimeMultTransSP)
	{
	    clReleaseProgram(program_multTransSP);
	    clReleaseKernel(kernel_multTransSP);
	}

	if (!isFirstTimeMultSP)
	{
		clReleaseMemObject(clDataSP);
		clReleaseMemObject(clLevelSP);
		clReleaseMemObject(clIndexSP);
	    clReleaseProgram(program_multSP);
	    clReleaseKernel(kernel_multSP);
	}

	isFirstTimeMultSP = true;
	isFirstTimeMultTransSP = true;
}

}
