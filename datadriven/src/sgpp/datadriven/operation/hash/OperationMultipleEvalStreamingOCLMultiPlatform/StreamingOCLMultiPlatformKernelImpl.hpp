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
#include <sgpp/base/opencl/OCLClonedBufferMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLConfigurationParameters.hpp>
#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLStretchedBufferMultiPlatform.hpp>
#include "StreamingOCLMultiPlatformKernelSourceBuilder.hpp"

namespace SGPP {
namespace datadriven {

template<typename T>
class StreamingOCLMultiPlatformKernelImpl {
private:
    size_t dims;

    cl_int err;

    base::OCLClonedBufferMultiPlatform deviceData;
    base::OCLClonedBufferMultiPlatform deviceLevel;
    base::OCLClonedBufferMultiPlatform deviceIndex;

    // use pinned memory (on host and device) to speed up data transfers from/to GPU
    base::OCLStretchedBufferMultiPlatform deviceGrid;
//    std::vector<T> hostGrid;
    base::OCLStretchedBufferMultiPlatform deviceTemp;
//    std::vector<T> hostTemp;

    std::map<cl_platform_id, std::vector<cl_kernel> > kernelsMultTrans;
    std::map<cl_platform_id, std::vector<cl_kernel> > kernelsMult;

    std::map<cl_platform_id, double *> deviceTimingsMult;
    std::map<cl_platform_id, double *> deviceTimingsMultTranspose;

    StreamingOCLMultiPlatformKernelSourceBuilder kernelSourceBuilder;
    std::shared_ptr<base::OCLManagerMultiPlatform> manager;
    std::shared_ptr<base::OCLConfigurationParameters> parameters;

    base::LinearLoadBalancerMultiPlatform multLoadBalancer;
    base::LinearLoadBalancerMultiPlatform multTransposeLoadBalancer;

    bool multKernelsBuilt = false;
    bool multTransKernelsBuilt = false;

public:

    StreamingOCLMultiPlatformKernelImpl(size_t dims, std::shared_ptr<base::OCLManagerMultiPlatform> manager,
            std::shared_ptr<base::OCLConfigurationParameters> parameters) :
            deviceData(manager), deviceLevel(manager), deviceIndex(manager), deviceGrid(manager), deviceTemp(manager), kernelSourceBuilder(
                    parameters, dims), manager(manager), parameters(parameters), multLoadBalancer(manager,
                    this->parameters), multTransposeLoadBalancer(manager, this->parameters) {

        this->dims = dims;

        for (base::OCLPlatformWrapper &platform : manager->platforms) {
            this->deviceTimingsMult[platform.platformId] = new double[platform.deviceCount];
            this->deviceTimingsMultTranspose[platform.platformId] = new double[platform.deviceCount];

            for (size_t i = 0; i < platform.deviceCount; i++) {
                //initialize with same timing to enforce equal problem sizes in the beginning
                this->deviceTimingsMult[platform.platformId][i] = 1.0;
                this->deviceTimingsMultTranspose[platform.platformId][i] = 1.0;
            }
        }

        this->err = CL_SUCCESS;

//        for (base::OCLPlatformWrapper &platform : manager->platforms) {
//            this->kernelsMult[platform.platformId] = new cl_kernel[platform.deviceCount];
//            this->kernelsMultTrans[platform.platformId] = new cl_kernel[platform.deviceCount];
//            // initialize arrays
//            for (size_t i = 0; i < platform.deviceCount; i++) {
//                this->kernelsMult[platform.platformId][i] = nullptr;
//                this->kernelsMultTrans[platform.platformId][i] = nullptr;
//            }
//        }
    }

    ~StreamingOCLMultiPlatformKernelImpl() {
        for (base::OCLPlatformWrapper &platform : manager->platforms) {
            if (kernelsMult.count(platform.platformId) > 0) {
                for (cl_kernel &kernelMult : kernelsMult[platform.platformId]) {
                    if (kernelMult != nullptr) {
                        clReleaseKernel(kernelMult);
                    }
                }
            }

            if (kernelsMultTrans.count(platform.platformId)) {
                for (cl_kernel &kernelMultTrans : kernelsMultTrans[platform.platformId]) {
                    if (kernelMultTrans != nullptr) {
                        clReleaseKernel(kernelMultTrans);
                    }
                }
            }

            delete[] this->deviceTimingsMult[platform.platformId];
            delete[] this->deviceTimingsMultTranspose[platform.platformId];
        }
    }

    void resetKernel() {
        //leads to a reallocation before next kernel execution
        releaseGridBuffers();
    }

    double mult(std::vector<T> &level, std::vector<T> &index, size_t gridSize, std::vector<T> &dataset, size_t datasetSize,
            std::vector<T> &alpha, std::vector<T> &result, const size_t start_index_grid, const size_t end_index_grid,
            const size_t start_index_data, const size_t end_index_data) {

        // check if there is something to do at all
        if (!(end_index_grid > start_index_grid && end_index_data > start_index_data)) {
            return 0.0;
        }

        double time = 0.0;

        if (!multKernelsBuilt) {
            std::string program_src = kernelSourceBuilder.generateSourceMult();
            manager->buildKernel(program_src, "multOCL", kernelsMult);
            multKernelsBuilt = true;
        }

        initOCLBuffers(level, index, gridSize, dataset, datasetSize);
        initParams(alpha, gridSize, result, datasetSize);

        // determine best fit
        std::map<cl_platform_id, size_t *> gpu_start_index_data;
        for (base::OCLPlatformWrapper &platform : manager->platforms) {
            gpu_start_index_data[platform.platformId] = new size_t[platform.deviceCount];
        }

        std::map<cl_platform_id, size_t *> gpu_end_index_data;
        for (base::OCLPlatformWrapper &platform : manager->platforms) {
            gpu_end_index_data[platform.platformId] = new size_t[platform.deviceCount];
        }

        multLoadBalancer.update(this->deviceTimingsMult);

        size_t dataBlockingSize = parameters->getAsUnsigned("KERNEL_DATA_BLOCKING_SIZE");
        multLoadBalancer.getPartitionSegments(start_index_data, end_index_data,
                parameters->getAsUnsigned("LOCAL_SIZE") * dataBlockingSize, gpu_start_index_data, gpu_end_index_data);

        // set kernel arguments
        cl_uint clResultSize = (cl_uint) datasetSize;
        cl_uint gpu_start_grid = (cl_uint) start_index_grid;
        cl_uint gpu_end_grid = (cl_uint) end_index_grid;
        //    std::cout << "start grid: " << gpu_start_grid << " end grid: " << gpu_end_grid << std::endl;

        for (base::OCLPlatformWrapper &platform : manager->platforms) {
            for (size_t i = 0; i < platform.deviceCount; i++) {
                cl_uint gpu_start_data = (cl_uint) gpu_start_index_data[platform.platformId][i]
                        / (cl_uint) dataBlockingSize;
                cl_uint gpu_end_data = (cl_uint) gpu_end_index_data[platform.platformId][i]
                        / (cl_uint) dataBlockingSize;

                cl_kernel kernel = kernelsMult[platform.platformId][i];

                if (gpu_end_data > gpu_start_data) {
                    err = clSetKernelArg(kernel, 0, sizeof(cl_mem),
                            this->deviceLevel.getBuffer(platform.platformId, i));
                    if (err != CL_SUCCESS) {
                        std::stringstream errorString;
                        errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                        throw SGPP::base::operation_exception(errorString.str());
                    }
                    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem),
                            this->deviceIndex.getBuffer(platform.platformId, i));
                    if (err != CL_SUCCESS) {
                        std::stringstream errorString;
                        errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                        throw SGPP::base::operation_exception(errorString.str());
                    }
                    err |= clSetKernelArg(kernel, 2, sizeof(cl_mem),
                            this->deviceData.getBuffer(platform.platformId, i));
                    if (err != CL_SUCCESS) {
                        std::stringstream errorString;
                        errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                        throw SGPP::base::operation_exception(errorString.str());
                    }
                    err |= clSetKernelArg(kernel, 3, sizeof(cl_mem),
                            this->deviceGrid.getBuffer(platform.platformId, i));
                    if (err != CL_SUCCESS) {
                        std::stringstream errorString;
                        errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                        throw SGPP::base::operation_exception(errorString.str());
                    }
                    err |= clSetKernelArg(kernel, 4, sizeof(cl_mem),
                            this->deviceTemp.getBuffer(platform.platformId, i));
                    if (err != CL_SUCCESS) {
                        std::stringstream errorString;
                        errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                        throw SGPP::base::operation_exception(errorString.str());
                    }
                    err |= clSetKernelArg(kernel, 5, sizeof(cl_uint), &clResultSize); // resultsize == number of entries in dataset
                    if (err != CL_SUCCESS) {
                        std::stringstream errorString;
                        errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                        throw SGPP::base::operation_exception(errorString.str());
                    }
                    err |= clSetKernelArg(kernel, 6, sizeof(cl_uint), &gpu_start_grid);
                    if (err != CL_SUCCESS) {
                        std::stringstream errorString;
                        errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                        throw SGPP::base::operation_exception(errorString.str());
                    }
                    err |= clSetKernelArg(kernel, 7, sizeof(cl_uint), &gpu_end_grid);
                    if (err != CL_SUCCESS) {
                        std::stringstream errorString;
                        errorString << "OCL Error: Failed to create kernel arguments for device " << i << std::endl;
                        throw SGPP::base::operation_exception(errorString.str());
                    }
                }
            }
        }

        std::map<cl_platform_id, cl_event *> clTimings;
        for (base::OCLPlatformWrapper &platform : manager->platforms) {
            clTimings[platform.platformId] = new cl_event[platform.deviceCount];
        }

        // enqueue kernel
        size_t local = parameters->getAsUnsigned("LOCAL_SIZE");

        std::map<cl_platform_id, size_t> platformTransferringDevice;

        bool allDistributed = false;
        for (base::OCLPlatformWrapper &platform : manager->platforms) {
            if (allDistributed) {
                platformTransferringDevice[platform.platformId] = 0;
                continue;
            }
            size_t activeDevices = 0;
            for (size_t i = 0; i < platform.deviceCount; i++) {
                size_t rangeSize = (gpu_end_index_data[platform.platformId][i] / dataBlockingSize)
                        - (gpu_start_index_data[platform.platformId][i] / dataBlockingSize);

                if (rangeSize > 0) {

                    size_t dataOffset = gpu_start_index_data[platform.platformId][i] / dataBlockingSize;
                    err = clEnqueueNDRangeKernel(platform.commandQueues[i], kernelsMult[platform.platformId][i], 1,
                            &dataOffset, &rangeSize, &local, 0, nullptr, &(clTimings[platform.platformId][i]));

                    if (err != CL_SUCCESS) {
                        std::stringstream errorString;
                        errorString << "OCL Error: Failed to enqueue kernel command! Error code: " << err << std::endl;
                        throw SGPP::base::operation_exception(errorString.str());
                    }

                    activeDevices++;
                } else {
                    allDistributed = true;
                    break;
                }
            }
            platformTransferringDevice[platform.platformId] = activeDevices;
        }

        //synchronization is realized through inorder queue in read

        deviceTemp.readFromBuffer(gpu_start_index_data, gpu_end_index_data);
        deviceTemp.combineBuffer(gpu_start_index_data, gpu_end_index_data, manager->platforms[0].platformId);

        T *hostTemp = (T *) deviceTemp.getMappedHostBuffer(manager->platforms[0].platformId);

        for (size_t i = start_index_data; i < end_index_data; i++) {
            result[i] = hostTemp[i];
        }

        // determine kernel execution time
        for (base::OCLPlatformWrapper &platform : manager->platforms) {
            for (size_t i = 0; i < platform.deviceCount; i++) {
                double tmpTime;
                cl_ulong startTime, endTime;
                startTime = endTime = 0;

                if (gpu_end_index_data[platform.platformId][i] > gpu_start_index_data[platform.platformId][i]) {
                    err = clGetEventProfilingInfo(clTimings[platform.platformId][i],
                    CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime, nullptr);

                    if (err != CL_SUCCESS) {
                        std::stringstream errorString;
                        errorString
                                << "OCL Error: Failed to read start-time from command queue (or crash in mult)! Error code: "
                                << err << std::endl;
                        throw SGPP::base::operation_exception(errorString.str());
                    }

                    err = clGetEventProfilingInfo(clTimings[platform.platformId][i],
                    CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, nullptr);

                    if (err != CL_SUCCESS) {
                        std::stringstream errorString;
                        errorString << "OCL Error: Failed to read end-time from command queue! Error code: " << err
                                << std::endl;
                        throw SGPP::base::operation_exception(errorString.str());
                    }
                }

                tmpTime = (double) (endTime - startTime);
                tmpTime *= 1e-9;
                this->deviceTimingsMult[platform.platformId][i] = tmpTime;

                if (tmpTime > time) {
                    time = tmpTime;
                }
            }

            // clean up
            for (size_t i = 0; i < platform.deviceCount; i++) {
                if (gpu_end_index_data[platform.platformId][i] > gpu_start_index_data[platform.platformId][i]) {
                    clReleaseEvent(clTimings[platform.platformId][i]);
                }
            }

            delete[] gpu_start_index_data[platform.platformId];
            delete[] gpu_end_index_data[platform.platformId];

            delete[] clTimings[platform.platformId];
        }

        return time;
    }

    double multTranspose(std::vector<T> &level, std::vector<T> &index, size_t gridSize, std::vector<T> &dataset,
            size_t datasetSize, std::vector<T> &source, std::vector<T> &result, const size_t start_index_grid,
            const size_t end_index_grid, const size_t start_index_data, const size_t end_index_data) {

        // check if there is something to do at all
        if (!(end_index_grid > start_index_grid && end_index_data > start_index_data)) {
            return 0.0;
        }

        double time = 0.0;
        size_t transGridBlockingSize = parameters->getAsUnsigned("KERNEL_TRANS_GRID_BLOCKING_SIZE");

        if (!multTransKernelsBuilt) {
            std::string program_src = kernelSourceBuilder.generateSourceMultTrans();
            manager->buildKernel(program_src, "multTransOCL", kernelsMultTrans);
            multTransKernelsBuilt = true;
        }

        initOCLBuffers(level, index, gridSize, dataset, datasetSize);
        initParams(result, gridSize, source, datasetSize);

        // determine best fit
        std::map<cl_platform_id, size_t *> gpu_start_index_grid;
        for (base::OCLPlatformWrapper &platform : manager->platforms) {
            gpu_start_index_grid[platform.platformId] = new size_t[platform.deviceCount];
        }

        std::map<cl_platform_id, size_t *> gpu_end_index_grid;
        for (base::OCLPlatformWrapper &platform : manager->platforms) {
            gpu_end_index_grid[platform.platformId] = new size_t[platform.deviceCount];
        }

        multTransposeLoadBalancer.update(this->deviceTimingsMultTranspose);
        multTransposeLoadBalancer.getPartitionSegments(start_index_grid, end_index_grid,
                parameters->getAsUnsigned("LOCAL_SIZE") * transGridBlockingSize, gpu_start_index_grid,
                gpu_end_index_grid);

        // set kernel arguments
        cl_uint clSourceSize = (cl_uint) datasetSize;
        cl_uint gpu_start_data = (cl_uint) start_index_data;
        cl_uint gpu_end_data = (cl_uint) end_index_data;

        for (base::OCLPlatformWrapper &platform : manager->platforms) {
            for (size_t i = 0; i < platform.deviceCount; i++) {
                cl_uint gpu_start_grid = (cl_uint) gpu_start_index_grid[platform.platformId][i]
                        / (cl_uint) transGridBlockingSize;
                cl_uint gpu_end_grid = (cl_uint) gpu_end_index_grid[platform.platformId][i]
                        / (cl_uint) transGridBlockingSize;

                cl_kernel kernel = kernelsMultTrans[platform.platformId][i];

                if (gpu_end_grid > gpu_start_grid) {
                    err = clSetKernelArg(kernel, 0, sizeof(cl_mem),
                            this->deviceLevel.getBuffer(platform.platformId, i));
                    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem),
                            this->deviceIndex.getBuffer(platform.platformId, i));
                    err |= clSetKernelArg(kernel, 2, sizeof(cl_mem),
                            this->deviceData.getBuffer(platform.platformId, i));
                    err |= clSetKernelArg(kernel, 3, sizeof(cl_mem),
                            this->deviceTemp.getBuffer(platform.platformId, i));
                    err |= clSetKernelArg(kernel, 4, sizeof(cl_mem),
                            this->deviceGrid.getBuffer(platform.platformId, i));
                    err |= clSetKernelArg(kernel, 5, sizeof(cl_uint), &clSourceSize); // sourceSize == number of entries in dataset
                    err |= clSetKernelArg(kernel, 6, sizeof(cl_uint), &gpu_start_data);
                    err |= clSetKernelArg(kernel, 7, sizeof(cl_uint), &gpu_end_data);
                    if (err != CL_SUCCESS) {
                        std::stringstream errorString;
                        errorString << "OCL Error: Failed to create kernel arguments for kernel " << i << std::endl;
                        throw SGPP::base::operation_exception(errorString.str());
                    }
                }
            }

        }

        std::map<cl_platform_id, cl_event *> clTimings;
        for (base::OCLPlatformWrapper &platform : manager->platforms) {
            clTimings[platform.platformId] = new cl_event[platform.deviceCount];
        }

        // enqueue kernels
        size_t local = parameters->getAsUnsigned("LOCAL_SIZE");

        std::map<cl_platform_id, size_t> platformTransferringDevice;

        bool allDistributed = false;
        for (base::OCLPlatformWrapper &platform : manager->platforms) {
            if (allDistributed) {
                platformTransferringDevice[platform.platformId] = 0;
                continue;
            }
            size_t activeDevices = 0;
            for (size_t i = 0; i < platform.deviceCount; i++) {
                size_t gpuEndGrid = gpu_end_index_grid[platform.platformId][i];
                size_t gpuStartGrid = gpu_start_index_grid[platform.platformId][i];
                size_t rangeSize = (gpuEndGrid / transGridBlockingSize) - (gpuStartGrid / transGridBlockingSize);

                if (rangeSize > 0) {

                    size_t gridOffset = gpu_start_index_grid[platform.platformId][i] / transGridBlockingSize;
                    err = clEnqueueNDRangeKernel(platform.commandQueues[i], kernelsMultTrans[platform.platformId][i], 1,
                            &gridOffset, &rangeSize, &local, 0, nullptr, &(clTimings[platform.platformId][i]));

                    if (err != CL_SUCCESS) {
                        std::stringstream errorString;
                        errorString << "OCL Error: Failed to enqueue kernel command! Error code: " << err << std::endl;
                        throw SGPP::base::operation_exception(errorString.str());
                    }

                    activeDevices++;
                } else {
                    allDistributed = true;
                    break;
                }
            }
            platformTransferringDevice[platform.platformId] = activeDevices;
        }

        //implicit synchonization in read due to inorder queue

        deviceGrid.readFromBuffer(gpu_start_index_grid, gpu_end_index_grid);
        deviceGrid.combineBuffer(gpu_start_index_grid, gpu_end_index_grid, manager->platforms[0].platformId);

        T *hostGrid = (T *) deviceGrid.getMappedHostBuffer(manager->platforms[0].platformId);

        for (size_t i = start_index_grid; i < end_index_grid; i++) {
//            printf("result i=%lu: %lf\n",i, hostGrid[i]);
            result[i] = hostGrid[i];
        }

        // determine kernel execution time
        for (base::OCLPlatformWrapper &platform : manager->platforms) {
            for (size_t i = 0; i < platform.deviceCount; i++) {
                double tmpTime;
                cl_ulong startTime, endTime;
                startTime = endTime = 0;

                if (gpu_end_index_grid[platform.platformId][i] > gpu_start_index_grid[platform.platformId][i]) {
                    err = clGetEventProfilingInfo(clTimings[platform.platformId][i],
                    CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime, nullptr);

                    if (err != CL_SUCCESS) {
                        std::stringstream errorString;
                        errorString
                                << "OCL Error: Failed to read start-time from command queue (or crash in multTranspose)! Error code: "
                                << err << std::endl;
                        throw SGPP::base::operation_exception(errorString.str());
                    }

                    err = clGetEventProfilingInfo(clTimings[platform.platformId][i],
                    CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, nullptr);

                    if (err != CL_SUCCESS) {
                        std::stringstream errorString;
                        errorString << "OCL Error: Failed to read end-time from command queue! Error code: " << err
                                << std::endl;
                        throw SGPP::base::operation_exception(errorString.str());
                    }
                }

                tmpTime = (double) (endTime - startTime);
                tmpTime *= 1e-9;
                this->deviceTimingsMultTranspose[platform.platformId][i] = tmpTime;

                if (tmpTime > time) {
                    time = tmpTime;
                }
            }

            // clean up
            for (size_t i = 0; i < platform.deviceCount; i++) {
                if (gpu_end_index_grid[platform.platformId][i] > gpu_start_index_grid[platform.platformId][i]) {
                    clReleaseEvent(clTimings[platform.platformId][i]);
                }
            }

            delete[] gpu_start_index_grid[platform.platformId];
            delete[] gpu_end_index_grid[platform.platformId];

            delete[] clTimings[platform.platformId];
        }

        return time;
    }
private:

    void releaseGridBuffers() {
        this->deviceLevel.freeBuffer();
        this->deviceIndex.freeBuffer();
        this->deviceGrid.freeBuffer();
    }

    void releaseDataBuffers() {
        this->deviceData.freeBuffer();
        this->deviceTemp.freeBuffer();
    }

    void initOCLBuffers(std::vector<T> level, std::vector<T> index, size_t gridSize, std::vector<T> dataset,
            size_t datasetSize) {
        if (!this->deviceLevel.isInitialized()) {
            this->deviceLevel.initializeBuffer(level.data(), sizeof(T), gridSize * this->dims);
        }

        if (!this->deviceIndex.isInitialized()) {
            this->deviceIndex.initializeBuffer(index.data(), sizeof(T), gridSize * this->dims);
        }

        if (!this->deviceData.isInitialized()) {
            this->deviceData.initializeBuffer(dataset.data(), sizeof(T), datasetSize * this->dims);
        }
    }

    void initParams(std::vector<T> grid, size_t gridSize, std::vector<T> tmp, size_t datasetSize) {
        if (!this->deviceGrid.isInitialized()) {
            this->deviceGrid.initializeBuffer(sizeof(T), gridSize * this->dims);
        }

        T *hostGrid = (T *) deviceGrid.getMappedHostBuffer(manager->platforms[0].platformId);
        for (size_t i = 0; i < gridSize; i++) {
            hostGrid[i] = grid[i];
        }

        deviceGrid.copyToOtherHostBuffers(manager->platforms[0].platformId);
        deviceGrid.writeToBuffer();

        if (!this->deviceTemp.isInitialized()) {
            this->deviceTemp.initializeBuffer(sizeof(T), datasetSize * this->dims);
        }

        T *hostTemp = (T*) this->deviceTemp.getMappedHostBuffer(manager->platforms[0].platformId);
        for (size_t i = 0; i < datasetSize; i++) {
            hostTemp[i] = tmp[i];
        }

        deviceTemp.copyToOtherHostBuffers(manager->platforms[0].platformId);
        deviceTemp.writeToBuffer();
    }

}
;

}
}
