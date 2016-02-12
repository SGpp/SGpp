// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <CL/cl.h>

#include <string.h>
#include <limits>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/LinearLoadBalancerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLClonedBufferSD.hpp>
#include "KernelSourceBuilderB.hpp"

namespace SGPP {
namespace datadriven {
namespace DensityOCLMultiPlatform {

template<typename T>
class KernelCreateGraph {
private:

	std::shared_ptr<base::OCLDevice> device;

	size_t dims;

	cl_int err;

	base::OCLClonedBufferSD<int> devicePoints;
	base::OCLClonedBufferSD<T> deviceData;
	base::OCLClonedBufferSD<T> deviceResultData;

	cl_kernel kernel;

	SGPP::datadriven::DensityOCLMultiPlatform::SourceBuilderCreateGraph<T> kernelSourceBuilder;

	std::shared_ptr<base::OCLManagerMultiPlatform> manager;

	double deviceTimingMult;

	json::Node &kernelConfiguration;


	bool verbose;

	size_t localSize;
	size_t dataBlockingSize;
	size_t scheduleSize;
	size_t totalBlockSize;
public:

	KernelCreateGraph(std::shared_ptr<base::OCLDevice> dev, size_t dims,
					  std::shared_ptr<base::OCLManagerMultiPlatform> manager, json::Node &kernelConfiguration) :
		device(dev), dims(dims), err(CL_SUCCESS), devicePoints(device),
		deviceData(device), deviceResultData(device), kernel(nullptr),
		kernelSourceBuilder(device, kernelConfiguration, dims), manager(manager), deviceTimingMult(0.0),
		kernelConfiguration(kernelConfiguration)
	{
		this->verbose = kernelConfiguration["VERBOSE"].getBool();

		if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0
		  && kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt() < dims) {
		  std::stringstream errorString;
		  errorString
		  << "OCL Error: setting \"KERNEL_DATA_STORE\" to \"register\" requires value of \"KERNEL_MAX_DIM_UNROLL\"";
		  errorString << " to be greater than the dimension of the data set, was set to"
		  << kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt() << "(device: \"" << device->deviceName
		  << "\")" << std::endl;
		  throw base::operation_exception(errorString.str());
		}

		localSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
		dataBlockingSize = kernelConfiguration["KERNEL_DATA_BLOCKING_SIZE"].getUInt();
		scheduleSize = kernelConfiguration["KERNEL_SCHEDULE_SIZE"].getUInt();
		totalBlockSize = dataBlockingSize * localSize;
	}

	~KernelCreateGraph()
	{
		if (kernel != nullptr) {
			clReleaseKernel(kernel);
			this->kernel = nullptr;
		}
	}

	void resetKernel()
	{
	}

	double create_graph(std::vector<T> &data, std::vector<T> &result, size_t k)
	{
		if (verbose)
		{
			std::cout << "entering mult, device: " << device->deviceName << " (" << device->deviceId << ")"
					  << std::endl;
		}

		//Build kernel if not already done
		if (this->kernel == nullptr)
		{
			if(verbose)
				std::cout<<"generating kernel source"<<std::endl;
			std::string program_src = kernelSourceBuilder.generateSource();
			if(verbose)
				std::cout<<"Source: "<<std::endl<<program_src<<std::endl;
			if(verbose)
				std::cout<<"building kernel"<<std::endl;
			this->kernel = manager->buildKernel(program_src, device, "cscheme");
		}

		//Load data into buffers if not already done
		if (!deviceDatao.isInitialized())
		{
			deviceData.intializeTo(data, 1, 0, data.size());
			std::vector<T> zeros(data.size()*k);
			for (size_t i = 0; i < data.size()*k; i++) {
				zeros[i] = 0.0;
			}
			deviceResultData.intializeTo(zeros, 1, 0, data.size()*k);
			clFinish(device->commandQueue);
		}
		this->deviceTimingMult = 0.0;
		size_t datasize=data.size()*k;

		//Set kernel arguments
		err = clSetKernelArg(this->kernel, 0, sizeof(cl_mem), this->deviceData.getBuffer());
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
			throw base::operation_exception(errorString.str());
		}
		err = clSetKernelArg(this->kernel, 1, sizeof(cl_uint), &dims);
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
			throw base::operation_exception(errorString.str());
		}
		err = clSetKernelArg(this->kernel, 2, sizeof(cl_uint), &k);
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
			throw base::operation_exception(errorString.str());
		}
		err = clSetKernelArg(this->kernel, 3, sizeof(cl_uint), &datasize);
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
			throw base::operation_exception(errorString.str());
		}
		err = clSetKernelArg(this->kernel, 4, sizeof(cl_mem), this->deviceResultData.getBuffer());
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
			throw base::operation_exception(errorString.str());
		}

		cl_event clTiming = nullptr;

		// enqueue kernel
		if(verbose)
			std::cout<<"Starting the kernel"<<std::endl;
		size_t *globalworkrange=new size_t[1];
		globalworkrange[0]=gridSize;
		err = clEnqueueNDRangeKernel(device->commandQueue, this->kernel, 1, 0, globalworkrange,
									 NULL, 0, nullptr, &clTiming);
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to enqueue kernel command! Error code: " << err << std::endl;
			throw base::operation_exception(errorString.str());
		}
		clFinish(device->commandQueue);

		if(verbose)
			std::cout<<"Finished kernel execution"<<std::endl;
		deviceResultData.readFromBuffer();
		clFinish(device->commandQueue);

		std::vector<T> &hostTemp = deviceResultData.getHostPointer();
		for(size_t i=0; i<gridSize; i++)
		{
			result[i]=hostTemp[i];
		}
		// determine kernel execution time
		cl_ulong startTime = 0;
		cl_ulong endTime = 0;

		err = clGetEventProfilingInfo(clTiming,
									  CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime, nullptr);

		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString
				<< "OCL Error: Failed to read start-time from command queue (or crash in mult)! Error code: "
				<< err << std::endl;
			throw base::operation_exception(errorString.str());
		}

		err = clGetEventProfilingInfo(clTiming,
									  CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, nullptr);

		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to read end-time from command queue! Error code: " << err
						<< std::endl;
			throw base::operation_exception(errorString.str());
		}

		clReleaseEvent(clTiming);

		double time = 0.0;
		time = (double) (endTime - startTime);
		time *= 1e-9;

		if (verbose)
		{
			{
				std::cout << "device: " << device->deviceName << " (" << device->deviceId << ") "
						  << "duration: " << time << std::endl;
			}
		}

		this->deviceTimingMult += time;
		return 0;
	}

};

}
}
}
