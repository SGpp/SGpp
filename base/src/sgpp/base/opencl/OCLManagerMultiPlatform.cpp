/*
 * OCLManager.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: pfandedd
 */
/*
#include <iostream>
#include <sstream>

#include "OCLManagerMultiPlatform.hpp"

#include <sgpp/base/exception/operation_exception.hpp>

namespace SGPP {
namespace base {

OCLManagerMultiPlatform::OCLManagerMultiPlatform(
		base::OpenCLConfigurationParameters parameters) :
		parameters(parameters), deviceType(0) {
	this->verbose = parameters.getAsBoolean("OCL_MANAGER_VERBOSE");

	cl_int err;

	this->setPlatformIDs();

	if (verbose) {
		this->printPlatformsInfo();
	}

	if (parameters["PLATFORM"].compare("all") == 0) {
		this->selectPlatform();
	}

	this->setDeviceType();

	this->setupDeviceIDs();

	if (verbose) {
		std::cout << "OCL Info: " << overallDeviceCount
				<< " OpenCL devices have been found!" << std::endl;
	}

	cl_uint maxDevices = (cl_uint) parameters.getAsUnsigned("MAX_DEVICES");

	if (maxDevices > this->overallDeviceCount) {
		std::stringstream errorString;
		errorString
				<< "OCL Error: More devices selected than available, currently can use at most "
				<< this->overallDeviceCount << " devices" << std::endl;
		throw SGPP::base::operation_exception(errorString.str());
	}

	if (parameters["SELECT_SPECIFIC_DEVICE"].compare("DISABLED") != 0) {
		if (maxDevices != 1) {
			std::stringstream errorString;
			errorString
					<< "OCL Error: Cannot select a specific device if more than one device is used, MAX_DEVICES be set incorrectly"
					<< std::endl;
			throw SGPP::base::operation_exception(errorString.str());
		}

		size_t selectedDevice = parameters.getAsUnsigned(
				"SELECT_SPECIFIC_DEVICE");

		if (selectedDevice > overallDeviceCount) {
			std::stringstream errorString;
			errorString
					<< "OCL Error: Illegal value set for \"SELECT_SPECIFIC_DEVICE\":"
					<< selectedDevice << std::endl;
			throw SGPP::base::operation_exception(errorString.str());
		}

		//always for the first platform
		//platformDeviceIds[0][0] = platformDeviceIds[0][selectedDevice];
		platforms[0].deviceIds[0] = platforms[0].deviceIds[selectedDevice];

		if (verbose) {
			std::cout << "OCL Info: select device number " << selectedDevice
					<< std::endl;
		}
	}

	if (maxDevices != 0 && maxDevices < overallDeviceCount) {
		overallDeviceCount = maxDevices;
	}

	if (verbose) {
		std::cout << "OCL Info: using " << overallDeviceCount << " device/s"
				<< std::endl;
	}

	for (size_t i = 0; i < platforms.size(); i++) {
//		cl_platform_id platformId = platformIds[i];
		cl_platform_id platformId = platforms[i].platformId;
		// Create OpenCL context
		cl_context_properties properties[3] = { CL_CONTEXT_PLATFORM,
				(cl_context_properties) platformId, 0 };

//		cl_context context = clCreateContext(properties, overallDeviceCount,
//				platformDeviceIds[platformId], nullptr, nullptr, &err);
		platforms[i].context = clCreateContext(properties, overallDeviceCount,
				platforms[i].deviceIds, nullptr, nullptr, &err);
//		this->platformContext.push_back(context);
	}

	if (err != CL_SUCCESS) {
		std::stringstream errorString;
		errorString
				<< "OCL Error: Failed to create OpenCL context! Error Code: "
				<< err << std::endl;
		throw SGPP::base::operation_exception(errorString.str());
	}

	// Creating the command queues
	for (size_t i = 0; i < platforms.size(); i++) {
		cl_platform_id platformId = platforms[i].platformId;
		for (size_t deviceIndex = 0; deviceIndex < platforms[i].deviceCount;
				deviceIndex++) {
			cl_device_id deviceId = platforms[i].deviceIds[deviceIndex];

			char buffer[128];
			err = clGetDeviceInfo(platforms[i].deviceIds[deviceIndex],
			CL_DEVICE_NAME, 128 * sizeof(char), &buffer, nullptr);

			if (err != CL_SUCCESS) {
				std::stringstream errorString;
				errorString
						<< "OCL Error: Failed to read the device name for device: "
						<< i << std::endl;
				throw SGPP::base::operation_exception(errorString.str());
			}

			if (verbose) {
				std::cout << "OCL Info: device name: " << buffer << std::endl;
			}

			platforms[i].commandQeues[deviceIndex] = clCreateCommandQueue(
					platforms[i].context, platforms[i].deviceIds[deviceIndex],
					CL_QUEUE_PROFILING_ENABLE, &err);

			if (err != CL_SUCCESS) {
				std::stringstream errorString;
				errorString
						<< "OCL Error: Failed to create command queue! Error Code: "
						<< err << std::endl;
				throw SGPP::base::operation_exception(errorString.str());
			}
		}

		if (verbose) {
			std::cout
					<< "OCL Info: Successfully initialized OpenCL (local workgroup size: "
					<< parameters.getAsUnsigned("LOCAL_SIZE") << ")"
					<< std::endl << std::endl;
		}
	}
}

/ **
 * @brief buildKernel builds the program that is represented by @a program_src and creates @a num_devices kernel objects
 * that are stored into the array @a kernel (must be already allocated with at least @a num_devices )
 *
 * @param program_src the source of the program to compile
 * @param kernel_name name of the kernel function (in program_src) to create the kernel for
 * @param context OpenCL context
 * @param num_devices number of OpenCL devices
 * @param device_ids array with device ids, necessary for displaying build info
 * @param kernel already allocated array: the resulting kernels are put into this array, one for each device (=> at least num_devices entries)
 * @return
 * /
void OCLManagerMultiPlatform::buildKernel(const std::string& program_src,
		const char* kernel_name, cl_context context, size_t num_devices,
		cl_device_id *device_ids, cl_kernel *kernel) {
	cl_int err;

	// setting the program
	const char* kernel_src = program_src.c_str();
	cl_program program = clCreateProgramWithSource(context, 1, &kernel_src,
	NULL, &err);

	if (err != CL_SUCCESS) {
		std::stringstream errorString;
		errorString << "OCL Error: Failed to create program! Error Code: "
				<< err << std::endl;
		throw SGPP::base::operation_exception(errorString.str());
	}

	std::string build_opts;

	if (parameters.getAsBoolean("ENABLE_OPTIMIZATIONS")) {
		//TODO: user should be able to change
		build_opts = parameters["OPTIMIZATION_FLAGS"];// -O5  -cl-mad-enable -cl-denorms-are-zero -cl-no-signed-zeros -cl-unsafe-math-optimizations -cl-finite-math-only -cl-fast-relaxed-math
	} else {
		build_opts = "-cl-opt-disable"; // -g
	}

	//TODO: check multi device support
	// compiling the program
	err = clBuildProgram(program, 0, NULL, build_opts.c_str(), NULL, NULL);

	if (err != CL_SUCCESS) {
		// get the build log
		size_t len;
		clGetProgramBuildInfo(program, device_ids[0], CL_PROGRAM_BUILD_LOG, 0,
		NULL, &len);
		std::string buffer(len, '\0');
		clGetProgramBuildInfo(program, device_ids[0], CL_PROGRAM_BUILD_LOG, len,
				&buffer[0], NULL);
		buffer = buffer.substr(0, buffer.find('\0'));

		if (verbose) {
			std::cout << "--- Build Log ---" << std::endl << buffer
					<< std::endl;
		}

		std::stringstream errorString;
		errorString << "OCL Error: OpenCL build error. Error code: " << err
				<< std::endl;
		throw SGPP::base::operation_exception(errorString.str());

	}

	// creating the kernel
	for (size_t i = 0; i < num_devices; i++) {
		kernel[i] = clCreateKernel(program, kernel_name, &err);

		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Failed to create kernel! Error code: "
					<< err << std::endl;
			throw SGPP::base::operation_exception(errorString.str());
		}
	}

	if (program) {
		clReleaseProgram(program);
	}
}

void OCLManagerMultiPlatform::setPlatformIDs() {
	cl_int err;

	// determine number of available OpenCL platforms
	cl_uint platformCount;
	err = clGetPlatformIDs(0, nullptr, &platformCount);

	if (err != CL_SUCCESS) {
		std::stringstream errorString;
		errorString
				<< "OCL Error: Unable to get number of OpenCL platforms. Error Code: "
				<< err << std::endl;
		throw SGPP::base::operation_exception(errorString.str());
	}

	if (verbose) {
		std::cout << "OCL Info: " << platforms.size()
				<< " OpenCL Platforms have been found" << std::endl;
	}

	// get available platforms
	cl_platform_id *platformIds = new cl_platform_id[platformCount];
	err = clGetPlatformIDs(platformCount, platformIds, nullptr);

	for (size_t i = 0; i < platforms.size(); i++) {
		OCLPlatformWrapper platformWrapper(platformIds[i]);
		platforms.push_back(platformWrapper);
	}

	if (err != CL_SUCCESS) {
		std::stringstream errorString;
		errorString << "OCL Error: Unable to get Platform ID. Error Code: "
				<< err << std::endl;
		throw SGPP::base::operation_exception(errorString.str());
	}
	delete platformIds;
}

void OCLManagerMultiPlatform::printPlatformsInfo() {
	cl_int err;

	for (cl_uint i = 0; i < platforms.size(); i++) {
		if (verbose) {
			char vendor_name[128] = { 0 };
			err = clGetPlatformInfo(platformIds[i], CL_PLATFORM_VENDOR,
					128 * sizeof(char), vendor_name, nullptr);

			if (CL_SUCCESS != err) {
				std::stringstream errorString;
				errorString << "OCL Error: Can't get platform vendor!"
						<< std::endl;
				throw SGPP::base::operation_exception(errorString.str());
			} else {
				std::cout << "OCL Info: Platform " << i << " vendor name: "
						<< vendor_name << std::endl;
			}
		}
	}
}

void OCLManagerMultiPlatform::selectPlatform() {
	cl_int err;

	for (cl_uint i = 0; i < platformCount; i++) {
		char platform_name[128] = { 0 };
		err = clGetPlatformInfo(platformIds[i], CL_PLATFORM_NAME,
				128 * sizeof(char), platform_name, nullptr);

		if (CL_SUCCESS != err) {
			std::stringstream errorString;
			errorString << "OCL Error: Can't get platform name!" << std::endl;
			throw SGPP::base::operation_exception(errorString.str());
		} else {
			if (verbose) {
				std::cout << "OCL Info: Platform " << i << " name: "
						<< platform_name << std::endl;
			}

			if (parameters["PLATFORM"].compare(platform_name) == 0) {
				platformIds[0] = platformIds[i];
				platformCount = 1;

				if (verbose) {
					std::cout << "platform selected" << std::endl;
				}
			}
		}
	}

	if (parameters["PLATFORM"].compare("first") == 0) {
		if (verbose) {
			std::cout << "using first platform" << std::endl;
		}

		platformCount = 1;
	}
}

void OCLManagerMultiPlatform::setTotalDeviceCount() {
	cl_int err;
	overallDeviceCount = 0;

	for (size_t i = 0; i < this->platformCount; i++) {
		cl_platform_id platformId = this->platformIds[i];
		uint32_t currentPlatformDevices;
		// get the number of devices
		err = clGetDeviceIDs(platformId, this->deviceType, 0, nullptr,
				&currentPlatformDevices);

		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Unable to get device count. Error Code: "
					<< err << std::endl;
			throw SGPP::base::operation_exception(errorString.str());
		}

		overallDeviceCount += currentPlatformDevices;
		this->platformDeviceCount[platformId] = currentPlatformDevices;
	}
}

void OCLManagerMultiPlatform::setupDeviceIDs() {
	cl_int err;
	size_t deviceIndex = 0;

	for (size_t i = 0; i < this->platformCount; i++) {
		cl_platform_id platformId = this->platformIds[i];
		cl_device_id *deviceIds =
				new cl_device_id[this->platformDeviceCount[platformId]];

		err = clGetDeviceIDs(platform_id, this->deviceType,
				this->platformDeviceCount[platformId], deviceIds, nullptr);

		this->platformDeviceIds[platformId] = deviceIds;

		if (err != CL_SUCCESS) {
			std::stringstream errorString;
			errorString << "OCL Error: Unable to get device id for platform "
					<< i << ". Error Code: " << err << std::endl;
			throw SGPP::base::operation_exception(errorString.str());
		}
	}
}

void OCLManagerMultiPlatform::setDeviceType() {
	if (parameters["DEVICE_TYPE"] == "CL_DEVICE_TYPE_CPU") {
		if (verbose) {
			std::cout << "OCL Info: looking for CPU device" << std::endl;
		}

		this->deviceType = CL_DEVICE_TYPE_CPU;
	} else if (parameters["DEVICE_TYPE"] == "CL_DEVICE_TYPE_GPU") {
		if (verbose) {
			std::cout << "OCL Info: looking for GPU device" << std::endl;
		}

		this->deviceType = CL_DEVICE_TYPE_GPU;
	} else if (parameters["DEVICE_TYPE"] == "CL_DEVICE_TYPE_ACCELERATOR") {
		if (verbose) {
			std::cout << "OCL Info: looking for device of accelerator type"
					<< std::endl;
		}

		this->deviceType = CL_DEVICE_TYPE_ACCELERATOR;
	} else if (parameters["DEVICE_TYPE"] == "CL_DEVICE_TYPE_ALL") {
		if (verbose) {
			std::cout << "OCL Info: looking for device of all available devices"
					<< std::endl;
		}

		this->deviceType = CL_DEVICE_TYPE_ALL;
	} else {
		throw SGPP::base::operation_exception(
				"OCL Error: No device found or unknown type specified (supported are: CPU, GPU, accelerator and all)");
	}
}

}
}
*/
