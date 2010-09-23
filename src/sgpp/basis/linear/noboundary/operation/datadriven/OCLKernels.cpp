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

// include OpenCL
#include "CL/cl.h"

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

	cl_int err;
	cl_platform_id platform_id;
	cl_device_id device_id;
	cl_uint num_platforms;
	cl_uint num_devices;
	cl_kernel kernel;
	cl_context context;
	cl_command_queue command_queue;
	cl_program program;
	cl_mem clAlpha, clData, clLevel, clIndex, clResult;

	std::stringstream stream_program_src;

	stream_program_src << "__kernel void multTransSPOCL(__global float* ptrAlpha," << std::endl;
	stream_program_src << "						__global float* ptrData," << std::endl;
	stream_program_src << "						__global float* ptrLevel," << std::endl;
	stream_program_src << "						__global float* ptrIndex," << std::endl;
	stream_program_src << "						__global float* ptrResult," << std::endl;
	stream_program_src << "						uint storageSize)" << std::endl;
	stream_program_src << "{" << std::endl;
	stream_program_src << "	int globalIdx = get_global_id(0);" << std::endl << std::endl;
	stream_program_src << "	float eval, index_calc, abs, last, localSupport;" << std::endl << std::endl;
	stream_program_src << "	float myResult = 0.0;" << std::endl << std::endl;
	stream_program_src << "	// Iterate over all grid points" << std::endl;
	stream_program_src << "	for(int j = 0; j < storageSize; j++)" << std::endl;
	stream_program_src << "	{" << std::endl;
	stream_program_src << "		float curSupport = ptrAlpha[j];" << std::endl << std::endl;
	for (size_t d = 0; d < dims; d++)
	{
		stream_program_src << "		eval = ((ptrLevel[(j*" << dims << ")+" << d << "]) * (ptrData[(" << d << "*" << result_size << ")+globalIdx]));" << std::endl;
		stream_program_src << "		index_calc = eval - (ptrIndex[(j*" << dims << ")+" << d << "]);" << std::endl;
		stream_program_src << "		abs = fabs(index_calc);" << std::endl;
		stream_program_src << "		last = 1.0 - abs;" << std::endl;
		stream_program_src << "		localSupport = fmax(last, 0.0);" << std::endl;
		stream_program_src << "		curSupport *= localSupport;" << std::endl << std::endl;
	}
	stream_program_src << "		myResult += curSupport;" << std::endl;
	stream_program_src << "	}" << std::endl;
	stream_program_src << "	ptrResult[globalIdx] = myResult;" << std::endl;
	stream_program_src << "}" << std::endl;

	std::string program_src = stream_program_src.str();

	//std::cout << program_src << std::endl;

	// Check for OCL Platforms
	err = clGetPlatformIDs(1, &platform_id, &num_platforms);
	if (err != CL_SUCCESS)
	{
		std::cout << "Unable to get Platform ID. Error Code: " << err << std::endl;
		return 0.0;
	}

	// Find out how many devices there are
	err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &num_devices);
    if (err != CL_SUCCESS)
    {
    	std::cout << "Unable to get Device ID. Error Code: " << err << std::endl;
    	return 0.0;
    }
    else if (num_devices == 0)
    {
    	std::cout << "NO GPU OpenCL devices have been found!" << std::endl;
    	return 0.0;
    }

    // Create GPU context
    context = clCreateContext(0, num_devices, &device_id, NULL, NULL, &err);
    if (err != CL_SUCCESS)
	{
    	std::cout << "Failed to create OpenCL context! Error Code: " << err << std::endl;
    	return 0.0;
	}

    // Creating the command queue
    command_queue = clCreateCommandQueue(context, device_id, 0, &err);
    if (err != CL_SUCCESS)
	{
    	std::cout << "Failed to create command queue! Error Code: " << err << std::endl;
    	return 0.0;
	}

    // setting the program
    const char* kernel_src = program_src.c_str();
    program = clCreateProgramWithSource(context, 1, &kernel_src, NULL, &err);
    if (err != CL_SUCCESS)
	{
    	std::cout << "Failed to create program! Error Code: " << err << std::endl;
    	return 0.0;
	}

    // compiling the program
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
    	std::cout << "OpenCL Build Error. Error Code: " << err << std::endl;

    	size_t len;
    	char buffer[2048];

    	// get the build log
    	clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);

    	std::cout << "--- Build Log ---" << std::endl << buffer << std::endl;
    	return 0.0;
    }

    // creating the kernel
    kernel = clCreateKernel(program, "multTransSPOCL", &err);
    if (err != CL_SUCCESS)
	{
    	std::cout << "Failed to create kernel! Error Code: " << err << std::endl;
    	return 0.0;
	}

    clAlpha = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*storageSize, ptrAlpha, NULL);
    clData = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*dims*result_size, ptrData, NULL);
    clLevel = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*dims*storageSize, ptrLevel, NULL);
    clIndex = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*dims*storageSize, ptrIndex, NULL);
    clResult = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float)*result_size, NULL, NULL);

    cl_uint clDims, clStorageSize, clResultSize;
    clDims = (cl_uint)dims;
    clStorageSize = (cl_uint)storageSize;
    clResultSize = (cl_uint)result_size;

    // set kernel arguments
    if ( clSetKernelArg(kernel, 0, sizeof(cl_mem), &clAlpha) ||
    		clSetKernelArg(kernel, 1, sizeof(cl_mem), &clData) ||
    		clSetKernelArg(kernel, 2, sizeof(cl_mem), &clLevel) ||
    		clSetKernelArg(kernel, 3, sizeof(cl_mem), &clIndex) ||
    		clSetKernelArg(kernel, 4, sizeof(cl_mem), &clResult) ||
    		clSetKernelArg(kernel, 5, sizeof(cl_uint), &clStorageSize) != CL_SUCCESS)
    {
    	std::cout << "Failed to create kernel Args!" << std::endl;
    	return 0.0;
    }

    size_t global[1];

    global[0] = result_size;

    // enqueue kernel
    err = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, global, NULL, 0, NULL, NULL);
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
    	std::cout << "Failed to enqueue read buffer command! Error Code: " << err << std::endl;
    	return 0.0;
	}

    // clean up
    clReleaseMemObject(clAlpha);
    clReleaseMemObject(clData);
    clReleaseMemObject(clLevel);
    clReleaseMemObject(clIndex);
    clReleaseMemObject(clResult);
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(command_queue);
    clReleaseContext(context);

//    for (size_t i = 0; i < result_size; i++)
//	{
//		for (size_t j = 0; j < storageSize; j++)
//		{
//			float curSupport = ptrAlpha[j];
//
//			for (size_t d = 0; d < dims; d++)
//			{
//				float eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*result_size)+i]));
//				float index_calc = eval - (ptrIndex[(j*dims)+d]);
//				float abs = fabs(index_calc);
//				float last = 1.0 - abs;
//				float localSupport = std::max<float>(last, 0.0);
//				curSupport *= localSupport;
//			}
//
//			ptrResult[i] += curSupport;
//		}
//	}

	return time;
}

}
