// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//#ifdef USEOCL
#include <omp.h>

#include <sgpp/globaldef.hpp>

#include "OCLKernelImpl.hpp"

namespace SGPP {
namespace datadriven {

void OCLKernelImpl::getPartitionSegment(size_t start, size_t end,
		size_t segmentCount, size_t segmentNumber, size_t* segmentStart,
		size_t* segmentEnd, size_t blockSize) {
	size_t totalSize = end - start;

	// check for valid input
	if (blockSize == 0) {
		throw SGPP::base::operation_exception("blockSize must not be zero!");
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

void OCLKernelImpl::getOpenMPPartitionSegment(size_t start,
		size_t end, size_t* segmentStart, size_t* segmentEnd,
		size_t blocksize) {
	size_t threadCount = omp_get_num_threads();
	size_t myThreadNum = omp_get_thread_num();
	getPartitionSegment(start, end, threadCount, myThreadNum, segmentStart,
			segmentEnd, blocksize);
}

OCLKernelImpl::OCLKernelImpl(OCLManager &manager): manager(manager) {

	this->num_devices = manager.num_devices;
	this->context = manager.context;
	this->command_queue = manager.command_queue;
	this->device_ids = manager.device_ids;
	this->OCLLocalSize = manager.getOCLLocalSize();

	this->clData = new cl_mem[num_devices];
	this->clLevel = new cl_mem[num_devices];
	this->clIndex = new cl_mem[num_devices];
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

OCLKernelImpl::~OCLKernelImpl() {
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

void OCLKernelImpl::releaseGridBuffers() {
	for (size_t i = 0; i < num_devices; i++) {
		if (clLevel[i]) {
			clReleaseMemObject(clLevel[i]);
			clLevel[i] = nullptr;
		}

		if (clIndex[i]) {
			clReleaseMemObject(clIndex[i]);
			clIndex[i] = nullptr;
		}

		if (clDevGrid[i]) {
			clReleaseMemObject(clDevGrid[i]);
			clDevGrid[i] = nullptr;
		}
	}

	if (clPinnedGrid) {
		clReleaseMemObject(clPinnedGrid);
		clPinnedGrid = nullptr;
	}
}

void OCLKernelImpl::releaseDataBuffers() {
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

void OCLKernelImpl::releaseKernelsAndPrograms() {
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

void OCLKernelImpl::resetKernel() {
	//leads to a reallocation before next kernel execution
	releaseGridBuffers();
}

}
}

//#endif //OPENCL
