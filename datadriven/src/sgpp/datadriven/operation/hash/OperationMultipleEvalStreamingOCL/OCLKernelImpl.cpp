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

unsigned int get_ocl_local_size() {
	// read environment variable for local work group size: SGPP_OCL_LOCAL_SIZE
	const char* ocl_local_size_env = getenv("SGPP_OCL_LOCAL_SIZE");

	if (ocl_local_size_env != NULL) {
		unsigned int num_ocl_devices_envvalue = (unsigned int) (strtoul(
				ocl_local_size_env, NULL, 0));

		if (num_ocl_devices_envvalue != 0) {
			return num_ocl_devices_envvalue;
		} else {
			std::cout << "Ignoring value: \"" << ocl_local_size_env
					<< "\" for SGPP_OCL_LOCAL_SIZE" << std::endl;
		}
	}

	return 64;
}

unsigned int OCLKernelImpl::ocl_local_size = get_ocl_local_size();

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

	clData = new cl_mem[num_devices];
	clLevel = new cl_mem[num_devices];
	clIndex = new cl_mem[num_devices];
	clDevGrid = new cl_mem[num_devices];
	clDevTmp = new cl_mem[num_devices];
	clPinnedGrid = NULL;
	clPinnedTmp = NULL;

	kernel_mult = new cl_kernel[num_devices];
	kernel_multTrans = new cl_kernel[num_devices];

	// initialize arrays
	for (size_t i = 0; i < num_devices; i++) {
		command_queue[i] = NULL;

		clData[i] = NULL;
		clLevel[i] = NULL;
		clIndex[i] = NULL;

		clDevGrid[i] = NULL;
		clDevTmp[i] = NULL;

		kernel_mult[i] = NULL;
		kernel_multTrans[i] = NULL;
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
			clLevel[i] = NULL;
		}

		if (clIndex[i]) {
			clReleaseMemObject(clIndex[i]);
			clIndex[i] = NULL;
		}

		if (clDevGrid[i]) {
			clReleaseMemObject(clDevGrid[i]);
			clDevGrid[i] = NULL;
		}
	}

	if (clPinnedGrid) {
		clReleaseMemObject(clPinnedGrid);
		clPinnedGrid = NULL;
	}
}

void OCLKernelImpl::releaseDataBuffers() {
	for (size_t i = 0; i < num_devices; i++) {
		if (clData[i]) {
			clReleaseMemObject(clData[i]);
			clData[i] = NULL;
		}

		if (clDevTmp[i]) {
			clReleaseMemObject(clDevTmp[i]);
			clDevTmp[i] = NULL;
		}
	}

	if (clPinnedTmp) {
		clReleaseMemObject(clPinnedTmp);
		clPinnedTmp = NULL;
	}
}

void OCLKernelImpl::releaseKernelsAndPrograms() {
	for (size_t i = 0; i < num_devices; i++) {
		if (kernel_mult[i]) {
			clReleaseKernel(kernel_mult[i]);
			kernel_mult[i] = NULL;
		}

		if (kernel_multTrans[i]) {
			clReleaseKernel(kernel_multTrans[i]);
			kernel_multTrans[i] = NULL;
		}
	}
}

void OCLKernelImpl::resetKernel() {
	releaseGridBuffers();
}

}
}

//#endif //OPENCL
