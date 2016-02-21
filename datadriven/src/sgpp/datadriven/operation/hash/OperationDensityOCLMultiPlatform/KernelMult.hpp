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
#include "KernelSourceBuilderMult.hpp"

namespace SGPP {
namespace datadriven {
namespace DensityOCLMultiPlatform {

template<typename T>
class KernelDensityMult {
private:

	std::shared_ptr<base::OCLDevice> device;

	size_t dims;
	size_t gridSize;
	T lambda;

	cl_int err;

	base::OCLClonedBufferSD<int> devicePoints;
	base::OCLClonedBufferSD<T> deviceAlpha;
	base::OCLClonedBufferSD<T> deviceResultData;

	cl_kernel kernelMult;

	SourceBuilderMult<T> kernelSourceBuilder;

	std::shared_ptr<base::OCLManagerMultiPlatform> manager;

	double deviceTimingMult;

	json::Node &kernelConfiguration;


	bool verbose;

	size_t localSize;
	size_t dataBlockingSize;
	size_t scheduleSize;
	size_t totalBlockSize;

public:

	KernelDensityMult(std::shared_ptr<base::OCLDevice> dev, size_t dims,
					  std::shared_ptr<base::OCLManagerMultiPlatform> manager, json::Node &kernelConfiguration,
					  std::vector<int> &points, T lambda) :
		device(dev), dims(dims), lambda(lambda), err(CL_SUCCESS), devicePoints(device),
		deviceAlpha(device), deviceResultData(device), kernelMult(nullptr),
		kernelSourceBuilder(device, kernelConfiguration, dims), manager(manager), deviceTimingMult(0.0),
		kernelConfiguration(kernelConfiguration)
		{
			this->verbose = true;
			gridSize = points.size()/(2*dims);
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

			devicePoints.intializeTo(points, 1, 0, gridSize*dims*2);
			clFinish(device->commandQueue);
		}

	~KernelDensityMult()
		{
			if (kernelMult != nullptr) {
				clReleaseKernel(kernelMult);
				this->kernelMult = nullptr;
			}
		}

	void resetKernel()
		{
		}

	double mult(std::vector<T> &alpha, std::vector<T> &result) {
		if (verbose)
		{
			std::cout << "entering mult, device: " << device->deviceName << " (" << device->deviceId << ")"
					  << std::endl;
		}

		//Build kernel if not already done
		if (this->kernelMult == nullptr)
		{
			if(verbose)
				std::cout<<"generating kernel source"<<std::endl;
			std::string program_src = kernelSourceBuilder.generateSource(dims, gridSize);
			if(verbose)
				std::cout<<"Source: "<<std::endl<<program_src<<std::endl;
			if(verbose)
				std::cout<<"building kernel"<<std::endl;
			this->kernelMult = manager->buildKernel(program_src, device, "multdensity");
		}

		//Load data into buffers if not already done
		deviceAlpha.intializeTo(alpha, 1, 0, gridSize);
		std::vector<T> zeros(gridSize);
		for (size_t i = 0; i < gridSize; i++) {
			zeros[i] = 0.0;
		}
		deviceResultData.intializeTo(zeros, 1, 0, gridSize);
		this->deviceTimingMult = 0.0;

		//Set kernel arguments
		err = clSetKernelArg(this->kernelMult, 0, sizeof(cl_mem), this->devicePoints.getBuffer());
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
			throw base::operation_exception(errorString.str());
		}
		err = clSetKernelArg(this->kernelMult, 1, sizeof(cl_mem), this->deviceAlpha.getBuffer());
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
			throw base::operation_exception(errorString.str());
		}
		err = clSetKernelArg(this->kernelMult, 2, sizeof(cl_mem), this->deviceResultData.getBuffer());
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
			throw base::operation_exception(errorString.str());
		}
		if(std::is_same<T, float>::value) {
			err = clSetKernelArg(this->kernelMult, 3, sizeof(cl_float), &lambda);
		}
		else {
			err = clSetKernelArg(this->kernelMult, 3, sizeof(cl_double), &lambda);
		}
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
		err = clEnqueueNDRangeKernel(device->commandQueue, this->kernelMult, 1, 0, globalworkrange,
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
	static void augmentDefaultParameters(SGPP::base::OCLOperationConfiguration &parameters) {
		for (std::string &platformName : parameters["PLATFORMS"].keys()) {
			json::Node &platformNode = parameters["PLATFORMS"][platformName];
			for (std::string &deviceName : platformNode["DEVICES"].keys()) {
				json::Node &deviceNode = platformNode["DEVICES"][deviceName];

				const std::string &kernelName = "multdensity";

				json::Node &kernelNode =
					deviceNode["KERNELS"].contains(kernelName) ?
					deviceNode["KERNELS"][kernelName] : deviceNode["KERNELS"].addDictAttr(kernelName);

				if (kernelNode.contains("VERBOSE") == false) {
					kernelNode.addIDAttr("VERBOSE", false);
				}

				if (kernelNode.contains("LOCAL_SIZE") == false) {
					kernelNode.addIDAttr("LOCAL_SIZE", 128ul);
				}

				if (kernelNode.contains("KERNEL_USE_LOCAL_MEMORY") == false) {
					kernelNode.addIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
				}

				if (kernelNode.contains("KERNEL_STORE_DATA") == false) {
					kernelNode.addTextAttr("KERNEL_STORE_DATA", "array");
				}

				if (kernelNode.contains("KERNEL_MAX_DIM_UNROLL") == false) {
					kernelNode.addIDAttr("KERNEL_MAX_DIM_UNROLL", 10ul);
				}

				if (kernelNode.contains("KERNEL_DATA_BLOCKING_SIZE") == false) {
					kernelNode.addIDAttr("KERNEL_DATA_BLOCKING_SIZE", 1ul);
				}

				if (kernelNode.contains("KERNEL_TRANS_GRID_BLOCKING_SIZE") == false) {
					kernelNode.addIDAttr("KERNEL_TRANS_GRID_BLOCKING_SIZE", 1ul);
				}

				if (kernelNode.contains("KERNEL_SCHEDULE_SIZE") == false) {
					kernelNode.addIDAttr("KERNEL_SCHEDULE_SIZE", 102400ul);
				}
			}
		}
	}
};

}
}
}
