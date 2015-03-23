// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <CL/cl.h>

#include <string.h>
#include <limits>

#include <sgpp/globaldef.hpp>

#include "OCLManager.hpp"
#include "OCLKernelSourceBuilder.hpp"

namespace SGPP {
namespace datadriven {

template<typename real_type>
class OCLKernelImpl {
private:
	size_t dims;

	cl_int err;
	cl_device_id* device_ids;
	cl_uint num_devices;
	cl_context context;
	cl_command_queue* command_queue;

	cl_mem* clData;
	cl_mem* clLevel;
	cl_mem* clIndex;

	// use pinned memory (on host and device) to speed up data transfers from/to GPU
	cl_mem* clDevGrid;
	cl_mem* clDevTmp;
	cl_mem clPinnedGrid;
	cl_mem clPinnedTmp;

	cl_kernel* kernel_multTrans;
	cl_kernel* kernel_mult;

	unsigned int OCLLocalSize;
	real_type* pinnedGrid;
	real_type* pinnedTmp;

	OCLKernelSourceBuilder<real_type> kernelSourceBuilder;
	OCLManager &manager;

public:

	OCLKernelImpl(size_t dims, OCLManager &manager) :
			manager(manager) {

		this->dims = dims;
		this->num_devices = manager.num_devices;
		this->context = manager.context;
		this->command_queue = manager.command_queue;
		this->device_ids = manager.device_ids;
		this->OCLLocalSize = manager.getOCLLocalSize();
		this->err = CL_SUCCESS;

		this->clData = new cl_mem[num_devices];
		this->clLevel = new cl_mem[num_devices];
		this->clIndex = new cl_mem[num_devices];
		this->pinnedGrid = nullptr;
		this->clDevGrid = new cl_mem[num_devices];
		this->clDevTmp = new cl_mem[num_devices];
		this->clPinnedGrid = nullptr;
		this->pinnedTmp = nullptr;
		this->clPinnedTmp = nullptr;

		this->kernel_mult = new cl_kernel[num_devices];
		this->kernel_multTrans = new cl_kernel[num_devices];

		// initialize arrays
		for (size_t i = 0; i < num_devices; i++) {
			this->clData[i] = nullptr;
			this->clLevel[i] = nullptr;
			this->clIndex[i] = nullptr;

			this->clDevGrid[i] = nullptr;
			this->clDevTmp[i] = nullptr;

			this->kernel_mult[i] = nullptr;
			this->kernel_multTrans[i] = nullptr;
		}
	}

	~OCLKernelImpl() {
		releaseDataBuffers();
		releaseGridBuffers();
		releaseKernelsAndPrograms();

		// release command queue
		for (size_t i = 0; i < num_devices; i++) {
			if (command_queue[i]) {
				clReleaseCommandQueue(command_queue[i]);
			}
		}

		// release context
		clReleaseContext(context);

		delete[] command_queue;

		delete[] clData;
		delete[] clLevel;
		delete[] clIndex;

		delete[] clDevGrid;
		delete[] clDevTmp;

		delete[] kernel_mult;
		delete[] kernel_multTrans;

		delete[] device_ids;
	}

	void resetKernel() {
		//leads to a reallocation before next kernel execution
		releaseGridBuffers();
	}

	double mult(real_type* level, real_type* index, size_t gridSize,
			real_type* dataset, size_t datasetSize, real_type* alpha,
			real_type* result, const size_t start_index_grid,
			const size_t end_index_grid, const size_t start_index_data,
			const size_t end_index_data) {

		// check if there is something to do at all
		if (!(end_index_grid > start_index_grid
				&& end_index_data > start_index_data)) {
			return 0.0;
		}

		double time = 0.0;

		if (kernel_mult[0] == nullptr) {
			this->createMult(this->dims, OCLLocalSize, context, num_devices,
					device_ids, kernel_mult);
		}

		initOCLBuffers(level, index, gridSize, dataset, datasetSize);
		initParams(alpha, gridSize, result, datasetSize);

		// determine best fit
		size_t* gpu_start_index_data = new size_t[num_devices];
		size_t* gpu_end_index_data = new size_t[num_devices];

		for (size_t gpu_num = 0; gpu_num < num_devices; gpu_num++) {
			this->getPartitionSegment(start_index_data, end_index_data,
					num_devices, gpu_num, &gpu_start_index_data[gpu_num],
					&gpu_end_index_data[gpu_num], OCLLocalSize);
		}

		// set kernel arguments
		cl_uint clResultSize = (cl_uint) datasetSize;
		cl_uint gpu_start_grid = (cl_uint) start_index_grid;
		cl_uint gpu_end_grid = (cl_uint) end_index_grid;

		for (size_t i = 0; i < num_devices; i++) {
			cl_uint gpu_start_data = (cl_uint) gpu_start_index_data[i];
			cl_uint gpu_end_data = (cl_uint) gpu_end_index_data[i];

			if (gpu_end_data > gpu_start_data) {
				if (clSetKernelArg(kernel_mult[i], 0, sizeof(cl_mem),
						&clLevel[i]) ||
						clSetKernelArg(kernel_mult[i], 1, sizeof(cl_mem), &clIndex[i]) ||
						clSetKernelArg(kernel_mult[i], 2, sizeof(cl_mem), &clData[i]) ||
						clSetKernelArg(kernel_mult[i], 3, sizeof(cl_mem), &clDevGrid[i]) ||
						clSetKernelArg(kernel_mult[i], 4, sizeof(cl_mem), &clDevTmp[i]) ||
						clSetKernelArg(kernel_mult[i], 5, sizeof(cl_uint), &clResultSize) || // resultsize == number of entries in dataset
						clSetKernelArg(kernel_mult[i], 6, sizeof(cl_uint), &gpu_start_grid) ||
						clSetKernelArg(kernel_mult[i], 7, sizeof(cl_uint), &gpu_end_grid) != CL_SUCCESS) {
					std::cout
							<< "OCL Error: Failed to create kernel Args for mult!"
							<< std::endl;
					return 0.0;
				}
			}
		}

		cl_event *clTimings = new cl_event[num_devices];
		cl_event *GPUDone = new cl_event[num_devices];

		// enqueue kernel
		size_t local = OCLLocalSize;
		size_t active_devices = 0;

		for (size_t i = 0; i < num_devices; i++) {
			size_t rangeSize = gpu_end_index_data[i] - gpu_start_index_data[i];

			if (rangeSize > 0) {
				err = clEnqueueNDRangeKernel(command_queue[i], kernel_mult[i],
						1, &gpu_start_index_data[i], &rangeSize, &local, 0,
						nullptr, &(clTimings[i]));

				if (active_devices != i) {
					std::cout
							<< "OCL Error: Splitting up calculations for multiple GPUs is erroneous, only the last chunks may handle 0 entries. syncing will be not correct.";
				}

				active_devices++;

				if (err != CL_SUCCESS) {
					std::cout
							<< "OCL Error: Failed to enqueue kernel command! Error Code: "
							<< err << std::endl;
					return 0.0;
				}
			}
		}

		// read data back
		for (size_t i = 0; i < num_devices; i++) {
			size_t rangeSize = gpu_end_index_data[i] - gpu_start_index_data[i];

			if (rangeSize > 0) {
				size_t offset = gpu_start_index_data[i];
				err = clEnqueueReadBuffer(command_queue[i], clDevTmp[i],
				CL_FALSE, sizeof(real_type) * offset,
						sizeof(real_type) * rangeSize, &(pinnedTmp[offset]), 0,
						nullptr, &(GPUDone[i]));

				if (err != CL_SUCCESS) {
					std::cout
							<< "OCL Error: Failed to enqueue read buffer command (mult)! Error Code: "
							<< err << std::endl;
					return 0.0;
				}
			}
		}

		// sync GPUs
		clWaitForEvents((cl_uint) active_devices, GPUDone);

		for (size_t i = start_index_data; i < end_index_data; i++) {
			result[i] = pinnedTmp[i];
		}

		// determine kernel execution time
		for (size_t i = 0; i < num_devices; i++) {
			double tmpTime;
			cl_ulong startTime, endTime;
			startTime = endTime = 0;

			if (gpu_end_index_data[i] > gpu_start_index_data[i]) {
				err = clGetEventProfilingInfo(clTimings[i],
				CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime,
						nullptr);

				if (err != CL_SUCCESS) {
					std::cout
							<< "OCL Error: Failed to read start-time from command queue! Error Code: "
							<< err << std::endl;
				}

				err = clGetEventProfilingInfo(clTimings[i],
				CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, nullptr);

				if (err != CL_SUCCESS) {
					std::cout
							<< "OCL Error: Failed to read end-time from command queue! Error Code: "
							<< err << std::endl;
				}
			}

			tmpTime = (double) (endTime - startTime);
			tmpTime *= 1e-9;

			if (tmpTime > time) {
				time = tmpTime;
			}
		}

		// clean up
		for (size_t i = 0; i < num_devices; i++) {
			if (gpu_end_index_data[i] > gpu_start_index_data[i]) {
				clReleaseEvent(clTimings[i]);
				clReleaseEvent(GPUDone[i]);
			}
		}

		delete[] gpu_start_index_data;
		delete[] gpu_end_index_data;

		delete[] clTimings;
		delete[] GPUDone;

		return time;
		//return 0.0;
	}

	double multTranspose(real_type* level, real_type* index, size_t gridSize,
			real_type* dataset, size_t datasetSize, real_type *source,
			real_type *result, const size_t start_index_grid,
			const size_t end_index_grid, const size_t start_index_data,
			const size_t end_index_data) {

		// check if there is something to do at all
		if (!(end_index_grid > start_index_grid
				&& end_index_data > start_index_data)) {
			return 0.0;
		}

		double time = 0.0;

		if (kernel_multTrans[0] == nullptr) {
			this->createMultTrans(this->dims, OCLLocalSize, context,
					num_devices, device_ids, kernel_multTrans);
		}

		initOCLBuffers(level, index, gridSize, dataset, datasetSize);
		initParams(result, gridSize, source, datasetSize);

		// determine best fit
		size_t* gpu_start_index_grid = new size_t[num_devices];
		size_t* gpu_end_index_grid = new size_t[num_devices];

		for (size_t gpu_num = 0; gpu_num < num_devices; gpu_num++) {
			this->getPartitionSegment(start_index_grid, end_index_grid,
					num_devices, gpu_num, &gpu_start_index_grid[gpu_num],
					&gpu_end_index_grid[gpu_num], OCLLocalSize);
		}

		// set kernel arguments
		cl_uint clSourceSize = (cl_uint) datasetSize;
		cl_uint gpu_start_data = (cl_uint) start_index_data;
		cl_uint gpu_end_data = (cl_uint) end_index_data;

		for (size_t i = 0; i < num_devices; i++) {
			cl_uint gpu_start_grid = (cl_uint) gpu_start_index_grid[i];
			cl_uint gpu_end_grid = (cl_uint) gpu_end_index_grid[i];

			if (gpu_end_grid > gpu_start_grid) {
				if (clSetKernelArg(kernel_multTrans[i], 0, sizeof(cl_mem),
						&clLevel[i]) ||
						clSetKernelArg(kernel_multTrans[i], 1, sizeof(cl_mem), &clIndex[i]) ||
						clSetKernelArg(kernel_multTrans[i], 2, sizeof(cl_mem), &clData[i]) ||
						clSetKernelArg(kernel_multTrans[i], 3, sizeof(cl_mem), &clDevTmp[i]) ||
						clSetKernelArg(kernel_multTrans[i], 4, sizeof(cl_mem), &clDevGrid[i]) ||
						clSetKernelArg(kernel_multTrans[i], 5, sizeof(cl_uint), &clSourceSize) || // sourceSize == number of entries in dataset
						clSetKernelArg(kernel_multTrans[i], 6, sizeof(cl_uint), &gpu_start_data) ||
						clSetKernelArg(kernel_multTrans[i], 7, sizeof(cl_uint), &gpu_end_data) != CL_SUCCESS) {
					std::cout
							<< "OCL Error: Failed to create kernel Args for kernel "
							<< i << "!" << std::endl;
					return 0.0;
				}
			}
		}

		cl_event* clTimings = new cl_event[num_devices];
		cl_event* GPUDone = new cl_event[num_devices];

		// enqueue kernels
		size_t local = OCLLocalSize;
		size_t active_devices = 0;

		for (size_t i = 0; i < num_devices; i++) {
			size_t rangeSize = gpu_end_index_grid[i] - gpu_start_index_grid[i];

			if (rangeSize > 0) {
				err = clEnqueueNDRangeKernel(command_queue[i],
						kernel_multTrans[i], 1, &gpu_start_index_grid[i],
						&rangeSize, &local, 0, nullptr, &(clTimings[i]));

				if (active_devices != i) {
					std::cout
							<< "OCL Error: Splitting up calculations for multiple GPUs is erroreous, only the last chunks may handle 0 entries. syncing will be not correct.";
				}

				active_devices++;

				if (err != CL_SUCCESS) {
					std::cout
							<< "OCL Error: Failed to enqueue kernel command! Error Code: "
							<< err << std::endl;
					return 0.0;
				}
			}
		}

		// read data back
		for (size_t i = 0; i < num_devices; i++) {
			size_t rangeSize = gpu_end_index_grid[i] - gpu_start_index_grid[i];

			if (rangeSize > 0) {
				size_t offset = gpu_start_index_grid[i];
				err = clEnqueueReadBuffer(command_queue[i], clDevGrid[i],
				CL_FALSE, sizeof(real_type) * offset,
						sizeof(real_type) * rangeSize, &(pinnedGrid[offset]), 0,
						nullptr, &(GPUDone[i]));

				if (err != CL_SUCCESS) {
					std::cout
							<< "OCL Error: Failed to enqueue read buffer command (multTrans)! Error Code: "
							<< err << std::endl;
					return 0.0;
				}
			}
		}

		// sync GPUs
		clWaitForEvents((cl_uint) active_devices, GPUDone);

		for (size_t i = start_index_grid; i < end_index_grid; i++) {
			result[i] = pinnedGrid[i];
		}

		// determine kernel execution time
		for (size_t i = 0; i < num_devices; i++) {
			double tmpTime;
			cl_ulong startTime, endTime;
			startTime = endTime = 0;

			if (gpu_end_index_grid[i] > gpu_start_index_grid[i]) {
				err = clGetEventProfilingInfo(clTimings[i],
				CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime,
						nullptr);

				if (err != CL_SUCCESS) {
					std::cout
							<< "OCL Error: Failed to read start-time from command queue! Error Code: "
							<< err << std::endl;
				}

				err = clGetEventProfilingInfo(clTimings[i],
				CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, nullptr);

				if (err != CL_SUCCESS) {
					std::cout
							<< "OCL Error: Failed to read end-time from command queue! Error Code: "
							<< err << std::endl;
				}
			}

			tmpTime = (double) (endTime - startTime);
			tmpTime *= 1e-9;

			if (tmpTime > time) {
				time = tmpTime;
			}
		}

		// clean up
		for (size_t i = 0; i < num_devices; i++) {
			if (gpu_end_index_grid[i] > gpu_start_index_grid[i]) {
				clReleaseEvent(clTimings[i]);
				clReleaseEvent(GPUDone[i]);
			}
		}

		delete[] gpu_start_index_grid;
		delete[] gpu_end_index_grid;

		delete[] clTimings;
		delete[] GPUDone;

		return time;
	}
private:

	cl_int createMultTrans(size_t dims, size_t local_workgroup_size,
			cl_context context, size_t num_devices, cl_device_id* device_ids,
			cl_kernel* kernel) {
		std::string program_src = kernelSourceBuilder.generateSourceMultTrans(
				dims, local_workgroup_size);
		return manager.buildKernel(program_src, "multTransOCL", context,
				num_devices, device_ids, kernel);
	}

	cl_int createMult(size_t dims, size_t local_workgroup_size,
			cl_context context, size_t num_devices, cl_device_id* device_ids,
			cl_kernel* kernel) {
		std::string program_src = kernelSourceBuilder.generateSourceMult(dims,
				local_workgroup_size);
		return manager.buildKernel(program_src, "multOCL", context, num_devices,
				device_ids, kernel);
	}

	size_t getChunkGridPoints() {
		return this->OCLLocalSize;
	}

	size_t getChunkDataPoints() {
		return this->OCLLocalSize;
	}

	void releaseGridBuffers() {
		for (size_t i = 0; i < num_devices; i++) {
			if (clLevel[i] != nullptr) {
				clReleaseMemObject(clLevel[i]);
				clLevel[i] = nullptr;
			}

			if (clIndex[i] != nullptr) {
				clReleaseMemObject(clIndex[i]);
				clIndex[i] = nullptr;
			}

			if (clDevGrid[i] != nullptr) {
				clReleaseMemObject(clDevGrid[i]);
				clDevGrid[i] = nullptr;
			}
		}

		if (clPinnedGrid) {
			clReleaseMemObject(clPinnedGrid);
			clPinnedGrid = nullptr;
		}
	}

	void releaseDataBuffers() {
		for (size_t i = 0; i < num_devices; i++) {
			if (clData[i]) {
				clReleaseMemObject(clData[i]);
				clData[i] = nullptr;
			}

			if (clDevTmp[i]) {
				clReleaseMemObject(clDevTmp[i]);
				clDevTmp[i] = nullptr;
			}
		}

		if (clPinnedTmp) {
			clReleaseMemObject(clPinnedTmp);
			clPinnedTmp = nullptr;
		}
	}

	void releaseKernelsAndPrograms() {
		for (size_t i = 0; i < num_devices; i++) {
			if (kernel_mult[i]) {
				clReleaseKernel(kernel_mult[i]);
				kernel_mult[i] = nullptr;
			}

			if (kernel_multTrans[i]) {
				clReleaseKernel(kernel_multTrans[i]);
				kernel_multTrans[i] = nullptr;
			}
		}
	}

	void initOCLBuffers(real_type* level, real_type* index, size_t gridSize,
			real_type* dataset, size_t datasetSize) {

		if (level != nullptr && clLevel[0] == nullptr) {
			for (size_t i = 0; i < num_devices; i++) {
				clLevel[i] = clCreateBuffer(context,
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
						sizeof(real_type) * gridSize * this->dims, level,
						nullptr);
			}
		}

		if (index != nullptr && clIndex[0] == nullptr) {
			for (size_t i = 0; i < num_devices; i++) {
				clIndex[i] = clCreateBuffer(context,
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
						sizeof(real_type) * gridSize * this->dims, index,
						nullptr);
			}
		}

		if (clData[0] == nullptr) { // use first element as indicator if data has been already copied to device
			for (size_t i = 0; i < num_devices; i++) {
				clData[i] = clCreateBuffer(context,
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
						sizeof(real_type) * datasetSize * this->dims, dataset,
						nullptr);
			}
		}
	}

	void initParams(real_type *grid, size_t gridSize, real_type * tmp,
			size_t datasetSize) {
		if (clPinnedGrid == nullptr) {
			size_t mem_size = sizeof(real_type) * gridSize * this->dims; //has to be the padded grid size
			clPinnedGrid = clCreateBuffer(context,
			CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, mem_size, nullptr,
					nullptr);

			for (size_t i = 0; i < num_devices; i++) {
				clDevGrid[i] = clCreateBuffer(context, CL_MEM_READ_WRITE,
						mem_size, nullptr, nullptr);
			}

			pinnedGrid = (real_type*) clEnqueueMapBuffer(command_queue[0],
					clPinnedGrid, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0,
					mem_size, 0, nullptr, nullptr, &err);
			if (err != CL_SUCCESS) {
				std::cout << "OCL Error: pinnedGrid: " << err << std::endl;
			}
		}

		for (size_t i = 0; i < gridSize; i++) {
			pinnedGrid[i] = grid[i];
		}

		for (size_t i = 0; i < num_devices; i++) {
			err = clEnqueueWriteBuffer(command_queue[i], clDevGrid[i], CL_TRUE,
					0, sizeof(real_type) * gridSize, pinnedGrid, 0, nullptr,
					nullptr);

			if (err != CL_SUCCESS) {
				std::cout
						<< "OCL Error: Failed to enqueue write buffer command (multTrans)! Error Code: "
						<< err << std::endl;
			}
		}

		if (clPinnedTmp == nullptr) {
			size_t mem_size = sizeof(real_type) * datasetSize * this->dims;
			clPinnedTmp = clCreateBuffer(context,
			CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, mem_size, nullptr,
					nullptr);

			for (size_t i = 0; i < num_devices; i++) {
				clDevTmp[i] = clCreateBuffer(context, CL_MEM_READ_WRITE,
						mem_size, nullptr, nullptr);
			}

			pinnedTmp = (real_type*) clEnqueueMapBuffer(command_queue[0],
					clPinnedTmp, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0,
					mem_size, 0, nullptr, nullptr, nullptr);
		}

		for (size_t i = 0; i < datasetSize; i++) {
			pinnedTmp[i] = tmp[i];
		}

		for (size_t i = 0; i < num_devices; i++) {
			clEnqueueWriteBuffer(command_queue[i], clDevTmp[i], CL_TRUE, 0,
					sizeof(real_type) * datasetSize, pinnedTmp, 0, nullptr,
					nullptr);
		}
	}

	void getPartitionSegment(size_t start, size_t end, size_t segmentCount,
			size_t segmentNumber, size_t* segmentStart, size_t* segmentEnd,
			size_t blockSize) {
		size_t totalSize = end - start;

		// check for valid input
		if (blockSize == 0) {
			throw SGPP::base::operation_exception(
					"blockSize must not be zero!");
		}

		if (totalSize % blockSize != 0) {
			//std::cout << "totalSize: " << totalSize << "; blockSize: " << blockSize << std::endl;
			throw SGPP::base::operation_exception(
					"totalSize must be divisible by blockSize without remainder, but it is not!");
		}

		// do all further calculations with complete blocks
		size_t blockCount = totalSize / blockSize;

		size_t blockSegmentSize = blockCount / segmentCount;
		size_t remainder = blockCount - blockSegmentSize * segmentCount;
		size_t blockSegmentOffset = 0;

		if (segmentNumber < remainder) {
			blockSegmentSize++;
			blockSegmentOffset = blockSegmentSize * segmentNumber;
		} else {
			blockSegmentOffset = remainder * (blockSegmentSize + 1)
					+ (segmentNumber - remainder) * blockSegmentSize;
		}

		*segmentStart = start + blockSegmentOffset * blockSize;
		*segmentEnd = *segmentStart + blockSegmentSize * blockSize;
	}

	void getOpenMPPartitionSegment(size_t start, size_t end,
			size_t* segmentStart, size_t* segmentEnd, size_t blocksize) {
		size_t threadCount = omp_get_num_threads();
		size_t myThreadNum = omp_get_thread_num();
		getPartitionSegment(start, end, threadCount, myThreadNum, segmentStart,
				segmentEnd, blocksize);
	}

};

}
}
