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
#include "KernelSourceBuilderCreateGraph.hpp"

namespace SGPP {
namespace datadriven {
namespace DensityOCLMultiPlatform {

template<typename T>
class KernelPruneGraph {
private:

	std::shared_ptr<base::OCLDevice> device;

	size_t dims;
	size_t k;

	cl_int err;

	base::OCLClonedBufferSD<int> devicePoints;
	base::OCLClonedBufferSD<T> deviceAlpha;
	base::OCLClonedBufferSD<T> deviceData;
	base::OCLClonedBufferSD<int> deviceGraph;

	cl_kernel kernel;

	SGPP::datadriven::DensityOCLMultiPlatform::SourceBuilderPruneGraph<T> kernelSourceBuilder;

	std::shared_ptr<base::OCLManagerMultiPlatform> manager;

	double deviceTimingMult;

	json::Node &kernelConfiguration;


	bool verbose;

	size_t localSize;
	size_t dataBlockingSize;
	size_t scheduleSize;
	size_t totalBlockSize;
public:

	KernelPruneGraph(std::shared_ptr<base::OCLDevice> dev, size_t dims, size_t k,
					  std::shared_ptr<base::OCLManagerMultiPlatform> manager, json::Node &kernelConfiguration) :
		device(dev), dims(dims), k(k), err(CL_SUCCESS),
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

	~KernelPruneGraph()
	{
		if (kernel != nullptr) {
			clReleaseKernel(kernel);
			this->kernel = nullptr;
		}
	}

	void resetKernel()
	{
	}

	double prune_graph(std::vector<int> &points, std::vector<T> &alpha, std::vector<T> &data, std::vector<int> &graph)
	{
		size_t gridSize = points.size()/(2*dims);
		size_t dataSize = data.size()/(i);
		if (verbose)
		{
			std::cout << "entering graph, device: " << device->deviceName << " (" << device->deviceId << ")"
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
			this->kernel = manager->buildKernel(program_src, device, "removeEdges");
		}

		//Load data into buffers if not already done
		if (!deviceData.isInitialized())
		{
			devicePoints.intializeTo(points, 1, 0, gridSize*dims*2);
			devicePoints.intializeTo(alpha, 1, 0, gridSize);
			deviceData.intializeTo(data, 1, 0, data.size());
			deviceGraph.intializeTo(graph, 1, 0, graph.size);
			clFinish(device->commandQueue);
		}
		this->deviceTimingMult = 0.0;
		size_t datasize=data.size()*k;

		//Set kernel arguments
		err = clSetKernelArg(this->kernel, 0, sizeof(cl_mem), this->deviceGraph.getBuffer());
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
			throw base::operation_exception(errorString.str());
		}
		err = clSetKernelArg(this->kernel, 1, sizeof(cl_mem), this->devicePoints.getBuffer());
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
			throw base::operation_exception(errorString.str());
		}
		err = clSetKernelArg(this->kernel, 2, sizeof(cl_mem), this->deviceData.getBuffer());
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
			throw base::operation_exception(errorString.str());
		}
		err = clSetKernelArg(this->kernel, 3, sizeof(cl_mem), this->deviceAlpha.getBuffer());
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
			throw base::operation_exception(errorString.str());
		}
		err = clSetKernelArg(this->kernel, 4, sizeof(cl_uint), &dims);
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
			throw base::operation_exception(errorString.str());
		}
		err = clSetKernelArg(this->kernel, 5, sizeof(cl_uint), &gridSize);
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
			throw base::operation_exception(errorString.str());
		}
		err = clSetKernelArg(this->kernel, 6, sizeof(cl_uint), &k);
		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
			throw base::operation_exception(errorString.str());
		}
		if(std::is_same<T, float>::value) {
		err = clSetKernelArg(this->kernelMult, 7, sizeof(cl_float), &lambda);
		}
		else {
		err = clSetKernelArg(this->kernelMult, 7, sizeof(cl_double), &lambda);
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
		globalworkrange[0]=data.size();
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

		std::vector<int> &hostTemp = deviceResultData.getHostPointer();
		for(size_t i=0; i<data.size()*k; i++)
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
				<< "OCL Error: Failed to read start-time from command queue (or crash in prune)! Error code: "
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
	static void augmentDefaultParameters(SGPP::base::OCLOperationConfiguration &parameters){
        for (std::string &platformName : parameters["PLATFORMS"].keys()) {
            json::Node &platformNode = parameters["PLATFORMS"][platformName];
            for (std::string &deviceName : platformNode["DEVICES"].keys()) {
                json::Node &deviceNode = platformNode["DEVICES"][deviceName];

                const std::string &kernelName = "connectNeighbors";

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
