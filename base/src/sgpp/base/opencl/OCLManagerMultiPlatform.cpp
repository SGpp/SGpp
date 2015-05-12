/*
 * OCLManager.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: pfandedd
 */

#include <iostream>
#include <sstream>

#include "OCLManagerMultiPlatform.hpp"

#include <sgpp/base/exception/operation_exception.hpp>

namespace SGPP {
namespace base {

OCLManagerMultiPlatform::OCLManagerMultiPlatform(base::OCLConfigurationParameters parameters) :
        parameters(parameters), deviceType(0) {
    this->verbose = parameters.getAsBoolean("OCL_MANAGER_VERBOSE");

    cl_int err;

    this->setPlatformIDs();

    if (verbose) {
        this->printPlatformsInfo();
    }

    this->setupPlatforms();

    this->setDeviceType();

    this->setTotalDeviceCount();

    this->setupDeviceIDs();

    if (verbose) {
        std::cout << "OCL Info: " << overallDeviceCount << " OpenCL devices have been found!" << std::endl;
    }

    if (parameters["SELECT_SPECIFIC_DEVICE"].compare("DISABLED") != 0) {
        if (this->platforms.size() > 1) {
            std::stringstream errorString;
            errorString << "OCL Error: Can only select a specific device if only one platform is used" << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }

        size_t selectedDevice = parameters.getAsUnsigned("SELECT_SPECIFIC_DEVICE");

        if (selectedDevice > platforms[0].deviceCount) {
            std::stringstream errorString;
            errorString << "OCL Error: Illegal value set for \"SELECT_SPECIFIC_DEVICE\":" << selectedDevice
                    << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }

        //always for the first platform
        platforms[0].deviceIds[0] = platforms[0].deviceIds[selectedDevice];
        platforms[0].deviceCount = 1;
        overallDeviceCount = 1;

        if (verbose) {
            std::cout << "OCL Info: select device number " << selectedDevice << std::endl;
        }

    }

//    if (maxDevices != 0 && maxDevices < overallDeviceCount) {
//        overallDeviceCount = maxDevices;
//    }

    if (verbose) {
        std::cout << "OCL Info: using " << overallDeviceCount << " device/s" << std::endl;
    }

    for (OCLPlatformWrapper &platform : platforms) {
        // Create OpenCL context
        cl_context_properties properties[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties) platform.platformId, 0 };

        platform.context = clCreateContext(properties, (cl_uint) platform.deviceCount, platform.deviceIds, nullptr,
                nullptr, &err);

        if (err != CL_SUCCESS) {
            std::stringstream errorString;
            errorString << "OCL Error: Failed to create OpenCL context! Error Code: " << err << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }
    }

    // Creating the command queues
    for (OCLPlatformWrapper &platform : platforms) {
        platform.commandQueues = new cl_command_queue[platform.deviceCount];
        for (size_t deviceIndex = 0; deviceIndex < platform.deviceCount; deviceIndex++) {

            char buffer[128];
            err = clGetDeviceInfo(platform.deviceIds[deviceIndex],
            CL_DEVICE_NAME, 128 * sizeof(char), &buffer, nullptr);

            if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to read the device name for device: " << deviceIndex << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
            }

            if (verbose) {
                std::cout << "OCL Info: device name: " << buffer << std::endl;
            }

            platform.commandQueues[deviceIndex] = clCreateCommandQueue(platform.context,
                    platform.deviceIds[deviceIndex],
                    CL_QUEUE_PROFILING_ENABLE, &err);

            if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create command queue! Error Code: " << err << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
            }
        }

        if (verbose) {
            std::cout << "OCL Info: Successfully initialized OpenCL (local workgroup size: "
                    << parameters.getAsUnsigned("LOCAL_SIZE") << ")" << std::endl << std::endl;
        }
    }
}

/**
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
 */
void OCLManagerMultiPlatform::buildKernel(const std::string &program_src, const char* kernel_name,
        std::map<cl_platform_id, cl_kernel *> &kernels) {
    cl_int err;

    for (OCLPlatformWrapper &platform : this->platforms) {
        // setting the program
        const char* kernel_src = program_src.c_str();
        cl_program program = clCreateProgramWithSource(platform.context, 1, &kernel_src,
        NULL, &err);

        if (err != CL_SUCCESS) {
            std::stringstream errorString;
            errorString << "OCL Error: Failed to create program! Error Code: " << err << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }

        std::string build_opts;

        if (parameters.getAsBoolean("ENABLE_OPTIMIZATIONS")) {
            //TODO: user should be able to change
            build_opts = parameters["OPTIMIZATION_FLAGS"]; // -O5  -cl-mad-enable -cl-denorms-are-zero -cl-no-signed-zeros -cl-unsafe-math-optimizations -cl-finite-math-only -cl-fast-relaxed-math
        } else {
            build_opts = "-cl-opt-disable"; // -g
        }

        //TODO: check multi device support
        // compiling the program
        err = clBuildProgram(program, 0, NULL, build_opts.c_str(), NULL, NULL);

        if (err != CL_SUCCESS) {
            // get the build log
            size_t len;
            clGetProgramBuildInfo(program, platform.deviceIds[0], CL_PROGRAM_BUILD_LOG, 0,
            NULL, &len);
            std::string buffer(len, '\0');
            clGetProgramBuildInfo(program, platform.deviceIds[0], CL_PROGRAM_BUILD_LOG, len, &buffer[0], NULL);
            buffer = buffer.substr(0, buffer.find('\0'));

            if (verbose) {
                std::cout << "--- Build Log ---" << std::endl << buffer << std::endl;
            }

            std::stringstream errorString;
            errorString << "OCL Error: OpenCL build error. Error code: " << err << std::endl;
            throw SGPP::base::operation_exception(errorString.str());

        }

        // creating the kernel
        for (size_t i = 0; i < platform.deviceCount; i++) {
            kernels[platform.platformId][i] = clCreateKernel(program, kernel_name, &err);

            if (err != CL_SUCCESS) {
                std::stringstream errorString;
                errorString << "OCL Error: Failed to create kernel! Error code: " << err << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
            }
        }

        if (program) {
            clReleaseProgram(program);
        }
    }
}

void OCLManagerMultiPlatform::setPlatformIDs() {
    cl_int err;

    // determine number of available OpenCL platforms
    cl_uint platformCount;
    err = clGetPlatformIDs(0, nullptr, &platformCount);

    if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Unable to get number of OpenCL platforms. Error Code: " << err << std::endl;
        throw SGPP::base::operation_exception(errorString.str());
    }

    if (verbose) {
        std::cout << "OCL Info: " << platforms.size() << " OpenCL Platforms have been found" << std::endl;
    }

    // get available platforms
    cl_platform_id *platformIds = new cl_platform_id[platformCount];
    err = clGetPlatformIDs(platformCount, platformIds, nullptr);

    for (size_t i = 0; i < platformCount; i++) {
        OCLPlatformWrapper platformWrapper(platformIds[i]);
        platforms.push_back(platformWrapper);
    }

    if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Unable to get Platform ID. Error Code: " << err << std::endl;
        throw SGPP::base::operation_exception(errorString.str());
    }
    delete platformIds;
}

void OCLManagerMultiPlatform::printPlatformsInfo() {
    cl_int err;

    for (cl_uint i = 0; i < platforms.size(); i++) {
        if (verbose) {
            char vendor_name[128] = { 0 };
            err = clGetPlatformInfo(platforms[i].platformId, CL_PLATFORM_VENDOR, 128 * sizeof(char), vendor_name,
                    nullptr);

            if (CL_SUCCESS != err) {
                std::stringstream errorString;
                errorString << "OCL Error: Can't get platform vendor!" << std::endl;
                throw SGPP::base::operation_exception(errorString.str());
            } else {
                std::cout << "OCL Info: Platform " << i << " vendor name: " << vendor_name << std::endl;
            }
        }
    }
}

void OCLManagerMultiPlatform::setupPlatforms() {
    cl_int err;

    size_t selectedPlatformIndex = 0;
    bool found = false;

    for (cl_uint i = 0; i < platforms.size(); i++) {
        char platform_name[128] = { 0 };
        err = clGetPlatformInfo(platforms[i].platformId, CL_PLATFORM_NAME, 128 * sizeof(char), platform_name, nullptr);

        if (CL_SUCCESS != err) {
            std::stringstream errorString;
            errorString << "OCL Error: Can't get platform name!" << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        } else {
            if (verbose) {
                std::cout << "OCL Info: Platform " << i << " name: " << platform_name << std::endl;
            }

            if (parameters["PLATFORM"].compare(platform_name) == 0) {
                selectedPlatformIndex = i;
                found = true;

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
        OCLPlatformWrapper selectedPlatform = platforms[0];
        platforms.clear();
        platforms.push_back(selectedPlatform);
    } else if (parameters["PLATFORM"].compare("all") != 0) {
        if (found) {
            OCLPlatformWrapper selectedPlatform = platforms[selectedPlatformIndex];
            platforms.clear();
            platforms.push_back(selectedPlatform);
        } else {
            std::stringstream errorString;
            errorString << "OCL Error: selected platform not found" << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }
    }

}

void OCLManagerMultiPlatform::setTotalDeviceCount() {
    cl_int err;
    overallDeviceCount = 0;

    for (size_t i = 0; i < platforms.size(); i++) {
        cl_platform_id platformId = platforms[i].platformId;
        uint32_t currentPlatformDevices;
        // get the number of devices
        err = clGetDeviceIDs(platformId, this->deviceType, 0, nullptr, &currentPlatformDevices);

        if (err != CL_SUCCESS) {
            std::stringstream errorString;
            errorString << "OCL Error: Unable to get device count. Error Code: " << err << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }

        overallDeviceCount += currentPlatformDevices;
        platforms[i].deviceCount = currentPlatformDevices;
    }
}

void OCLManagerMultiPlatform::setupDeviceIDs() {
    cl_int err;

    for (size_t i = 0; i < platforms.size(); i++) {
        cl_platform_id platformId = platforms[i].platformId;
        cl_device_id *deviceIds = new cl_device_id[platforms[i].deviceCount];

        err = clGetDeviceIDs(platformId, this->deviceType, (cl_uint) platforms[i].deviceCount, deviceIds, nullptr);

        platforms[i].deviceIds = deviceIds;

        if (err != CL_SUCCESS) {
            std::stringstream errorString;
            errorString << "OCL Error: Unable to get device id for platform " << i << ". Error Code: " << err
                    << std::endl;
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
            std::cout << "OCL Info: looking for device of accelerator type" << std::endl;
        }

        this->deviceType = CL_DEVICE_TYPE_ACCELERATOR;
    } else if (parameters["DEVICE_TYPE"] == "CL_DEVICE_TYPE_ALL") {
        if (verbose) {
            std::cout << "OCL Info: looking for device of all available devices" << std::endl;
        }

        this->deviceType = CL_DEVICE_TYPE_ALL;
    } else {
        throw SGPP::base::operation_exception(
                "OCL Error: No device found or unknown type specified (supported are: CPU, GPU, accelerator and all)");
    }
}

}
}
