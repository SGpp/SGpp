// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include "StreamingModOCLMaskMultiPlatformKernelImpl.hpp"

#include <string.h>
#include <limits>
#include <chrono>
#include <CL/cl.h>
#include <omp.h>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/LinearLoadBalancerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLClonedBufferSD.hpp>
#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLStretchedBuffer.hpp>
#include "StreamingModOCLMaskMultiPlatformKernelSourceBuilder.hpp"

namespace SGPP {
namespace datadriven {

template<typename real_type>
class StreamingModOCLMaskMultiPlatformKernelImpl {
private:

    std::shared_ptr<base::OCLDevice> device;

    size_t dims;

    cl_int err;

    base::OCLClonedBufferSD<real_type> deviceLevel;
    base::OCLClonedBufferSD<real_type> deviceIndex;
    base::OCLClonedBufferSD<real_type> deviceMask;
    base::OCLClonedBufferSD<real_type> deviceOffset;
    base::OCLClonedBufferSD<real_type> deviceAlpha;

    base::OCLClonedBufferSD<real_type> deviceData;
    base::OCLClonedBufferSD<real_type> deviceResultData;
    base::OCLClonedBufferSD<real_type> deviceResultGrid;

    cl_kernel kernelMult;
    cl_kernel kernelMultTranpose;

    double deviceTimingMult;
    double deviceTimingMultTranspose;

    StreamingModOCLMaskMultiPlatformKernelSourceBuilder kernelSourceBuilder;
    std::shared_ptr<base::OCLManagerMultiPlatform> manager;
    std::shared_ptr<base::OCLOperationConfiguration> parameters;

    std::shared_ptr<base::QueueLoadBalancer> queueLoadBalancerMult;
    std::shared_ptr<base::QueueLoadBalancer> queueLoadBalancerMultTranspose;

public:

    StreamingModOCLMaskMultiPlatformKernelImpl(std::shared_ptr<base::OCLDevice> device, size_t dims,
            std::shared_ptr<base::OCLManagerMultiPlatform> manager,
            std::shared_ptr<base::OCLOperationConfiguration> parameters,
            std::shared_ptr<base::QueueLoadBalancer> queueBalancerMult,
            std::shared_ptr<base::QueueLoadBalancer> queueBalancerMultTranpose) :
            device(device), dims(dims), err(CL_SUCCESS), deviceLevel(device), deviceIndex(device), deviceMask(device), deviceOffset(
                    device), deviceAlpha(device), deviceData(device), deviceResultData(device), deviceResultGrid(
                    device), kernelMult(nullptr), kernelMultTranpose(nullptr), kernelSourceBuilder(parameters, dims), manager(
                    manager), parameters(parameters), queueLoadBalancerMult(queueBalancerMult), queueLoadBalancerMultTranspose(
                    queueBalancerMultTranpose) {

        //initialize with same timing to enforce equal problem sizes in the beginning
        this->deviceTimingMult = 1.0;
        this->deviceTimingMultTranspose = 1.0;
    }

    ~StreamingModOCLMaskMultiPlatformKernelImpl() {
        if (this->kernelMult != nullptr) {
            clReleaseKernel(this->kernelMult);
            this->kernelMult = nullptr;
        }

        if (this->kernelMultTranpose != nullptr) {
            clReleaseKernel(this->kernelMultTranpose);
            this->kernelMultTranpose = nullptr;
        }
    }

    void resetKernel() {
        //leads to a reallocation before next kernel execution
        releaseGridBuffers();
    }

    double mult(std::vector<real_type> &level, std::vector<real_type> &index, std::vector<real_type> &mask,
            std::vector<real_type> &offset, size_t gridSize, std::vector<real_type> &dataset, size_t datasetSize,
            std::vector<real_type> &alpha, std::vector<real_type> &result, const size_t start_index_grid,
            const size_t end_index_grid, const size_t start_index_data, const size_t end_index_data) {

        // check if there is something to do at all
        if (!(end_index_grid > start_index_grid && end_index_data > start_index_data)) {
            return 0.0;
        }

        if (this->kernelMult == nullptr) {
            std::string program_src = kernelSourceBuilder.generateSourceMult();
            this->kernelMult = manager->buildKernel(program_src, device, "multOCLMask");

        }

        this->deviceTimingMult = 0.0;

        while (true) {

            size_t kernelStartData;
            size_t kernelEndData;

            // set kernel arguments
            size_t kernelStartGrid = start_index_grid;
            size_t kernelEndGrid = end_index_grid;

            //TODO: change after blocking is implemented
            //TODO: don't forget to set padding to DATA_BLOCKING * THREAD_BLOCK_SIZE
            size_t dataBlockingSize = 1;

            //TODO: start_index_data not considered!
            bool segmentAvailable = queueLoadBalancerMult->getNextSegment(kernelStartData, kernelEndData);
            if (!segmentAvailable) {
                break;
            }

            size_t rangeSizeUnblocked = kernelEndData - kernelStartData;

            std::cout << "device: " << device->deviceId << "kernel from: " << kernelStartData << " to: "
                    << kernelEndData << " -> range: " << rangeSizeUnblocked << std::endl;

            initGridBuffers(level, index, mask, offset, alpha, kernelStartGrid, kernelEndGrid);
            initDatasetBuffers(dataset, kernelStartData, kernelEndData);
            initDatasetResultBuffers(kernelStartData, kernelEndData);

            clFinish(device->commandQueue);
            std::cout << "wrote to device: " << device->deviceId << "" << std::endl;

            size_t rangeSizeBlocked = (kernelEndData / dataBlockingSize) - (kernelStartData / dataBlockingSize);

            if (rangeSizeBlocked > 0) {

                err = clSetKernelArg(this->kernelMult, 0, sizeof(cl_mem), this->deviceLevel.getBuffer());
                if (err != CL_SUCCESS) {
                    std::stringstream errorString;
                    errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
                    throw SGPP::base::operation_exception(errorString.str());
                }
                err = clSetKernelArg(this->kernelMult, 1, sizeof(cl_mem), this->deviceIndex.getBuffer());
                if (err != CL_SUCCESS) {
                    std::stringstream errorString;
                    errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
                    throw SGPP::base::operation_exception(errorString.str());
                }
                err = clSetKernelArg(this->kernelMult, 2, sizeof(cl_mem), this->deviceMask.getBuffer());
                if (err != CL_SUCCESS) {
                    std::stringstream errorString;
                    errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
                    throw SGPP::base::operation_exception(errorString.str());
                }
                err = clSetKernelArg(this->kernelMult, 3, sizeof(cl_mem), this->deviceOffset.getBuffer());
                if (err != CL_SUCCESS) {
                    std::stringstream errorString;
                    errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
                    throw SGPP::base::operation_exception(errorString.str());
                }
                err = clSetKernelArg(this->kernelMult, 4, sizeof(cl_mem), this->deviceData.getBuffer());
                if (err != CL_SUCCESS) {
                    std::stringstream errorString;
                    errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
                    throw SGPP::base::operation_exception(errorString.str());
                }
                err = clSetKernelArg(this->kernelMult, 5, sizeof(cl_mem), this->deviceAlpha.getBuffer());
                if (err != CL_SUCCESS) {
                    std::stringstream errorString;
                    errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
                    throw SGPP::base::operation_exception(errorString.str());
                }
                err = clSetKernelArg(this->kernelMult, 6, sizeof(cl_mem), this->deviceResultData.getBuffer());
                if (err != CL_SUCCESS) {
                    std::stringstream errorString;
                    errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
                    throw SGPP::base::operation_exception(errorString.str());
                }
                err = clSetKernelArg(this->kernelMult, 7, sizeof(cl_uint), &rangeSizeUnblocked); // resultsize == number of entries in dataset
                if (err != CL_SUCCESS) {
                    std::stringstream errorString;
                    errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
                    throw SGPP::base::operation_exception(errorString.str());
                }
                err = clSetKernelArg(this->kernelMult, 8, sizeof(cl_uint), &kernelStartGrid);
                if (err != CL_SUCCESS) {
                    std::stringstream errorString;
                    errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
                    throw SGPP::base::operation_exception(errorString.str());
                }
                err = clSetKernelArg(this->kernelMult, 9, sizeof(cl_uint), &kernelEndGrid);
                if (err != CL_SUCCESS) {
                    std::stringstream errorString;
                    errorString << "OCL Error: Failed to create kernel arguments for device " << std::endl;
                    throw SGPP::base::operation_exception(errorString.str());
                }

                cl_event clTiming = nullptr;

                // enqueue kernel
                size_t localSize = (*parameters)["LOCAL_SIZE"].getUInt();

                std::cout << "commandQueue: " << device->commandQueue << std::endl;

                char deviceName[128] = { 0 };
                cl_uint err;
                err = clGetDeviceInfo(device->deviceId,
                CL_DEVICE_NAME, 128 * sizeof(char), &deviceName, nullptr);

                if (err != CL_SUCCESS) {
                    std::stringstream errorString;
                    errorString << "OCL Error: Failed to read the device name for device: " << device->deviceId << std::endl;
                    throw SGPP::base::operation_exception(errorString.str());
                }

                std::cout << "OCL Info: detected device, name: \"" << deviceName << "\"" << std::endl;


                err = clEnqueueNDRangeKernel(device->commandQueue, this->kernelMult, 1, 0, &rangeSizeUnblocked,
                        &localSize, 0, nullptr, &clTiming);

                if (err != CL_SUCCESS) {
                    std::stringstream errorString;
                    errorString << "OCL Error: Failed to enqueue kernel command! Error code: " << err << std::endl;
                    throw SGPP::base::operation_exception(errorString.str());
                }

                clFinish(device->commandQueue);
                std::cout << "executed kernel: " << device->deviceId << "" << std::endl;

                //TODO: implement treatment of start_index_data in queueLoadBalancer!
                deviceResultData.readFromBuffer();
                std::cout << "read from device: " << device->deviceId << "" << std::endl;

                clFinish(device->commandQueue);

                std::vector<real_type> &hostTemp = deviceResultData.getHostPointer();
                size_t deviceIndex = 0;
                for (size_t i = 0; i < rangeSizeUnblocked; i++) {
//                    std::cout << "resultDevice[" << deviceIndex << "] = " << hostTemp[deviceIndex] << std::endl;
                    result[kernelStartData + i] = hostTemp[deviceIndex];
                    deviceIndex += 1;
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
                    throw SGPP::base::operation_exception(errorString.str());
                }

                err = clGetEventProfilingInfo(clTiming,
                CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, nullptr);

                if (err != CL_SUCCESS) {
                    std::stringstream errorString;
                    errorString << "OCL Error: Failed to read end-time from command queue! Error code: " << err
                            << std::endl;
                    throw SGPP::base::operation_exception(errorString.str());
                }

                clReleaseEvent(clTiming);

                double time = 0.0;
                time = (double) (endTime - startTime);
                time *= 1e-9;
                this->deviceTimingMult += time;
            }

        }

        return this->deviceTimingMult;
    }

    double multTranspose(std::vector<real_type> &level, std::vector<real_type> &index, std::vector<real_type> &mask,
            std::vector<real_type> &offset, size_t gridSize, std::vector<real_type> &dataset, size_t datasetSize,
            std::vector<real_type> &source, std::vector<real_type> &result, const size_t start_index_grid,
            const size_t end_index_grid, const size_t start_index_data, const size_t end_index_data) {
        /*
         // check if there is something to do at all
         if (!(end_index_grid > start_index_grid && end_index_data > start_index_data)) {
         return 0.0;
         }

         double time = 0.0;

         if (this->kernelMultTranpose == nullptr) {
         std::string program_src = kernelSourceBuilder.generateSourceMultTrans();
         manager->buildKernel(program_src, this->context, "multOCLMask",
         this->platformId, this->deviceId, this->kernelMultTranpose);
         }

         initOCLBuffers(level, index, mask, offset, gridSize, dataset, datasetSize);
         initParams(result, gridSize, source, datasetSize);

         // determine best fit
         size_t* gpu_start_index_grid = new size_t[manager->num_devices];
         size_t* gpu_end_index_grid = new size_t[manager->num_devices];

         multTransposeLoadBalancer.update(this->deviceTimingsMultTranspose);
         multTransposeLoadBalancer.getPartitionSegments(start_index_grid, end_index_grid,
         (*parameters)["LOCAL_SIZE"].getUInt(), gpu_start_index_grid,
         gpu_end_index_grid);

         // set kernel arguments
         cl_uint clSourceSize = (cl_uint) datasetSize;
         cl_uint gpu_start_data = (cl_uint) start_index_data;
         cl_uint gpu_end_data = (cl_uint) end_index_data;

         //    std::cout << "start data: " << gpu_start_data << " end data: " << gpu_end_data << std::endl;

         for (size_t i = 0; i < manager->num_devices; i++) {
         cl_uint gpu_start_grid = (cl_uint) gpu_start_index_grid[i];
         cl_uint gpu_end_grid = (cl_uint) gpu_end_index_grid[i];
         //      std::cout << "device: " << i << " start grid: " << gpu_start_grid << " end grid: " << gpu_end_grid << std::endl;

         if (gpu_end_grid > gpu_start_grid) {
         err = clSetKernelArg(kernel_multTrans[i], 0, sizeof(cl_mem),
         this->deviceLevel.getBuffer(i));
         if (err != CL_SUCCESS) {
         std::stringstream errorString;
         errorString
         << "OCL Error: Failed to create kernel arguments for device "
         << i << std::endl;
         throw SGPP::base::operation_exception(errorString.str());
         }
         err = clSetKernelArg(kernel_multTrans[i], 1, sizeof(cl_mem),
         this->deviceIndex.getBuffer(i));
         if (err != CL_SUCCESS) {
         std::stringstream errorString;
         errorString
         << "OCL Error: Failed to create kernel arguments for device "
         << i << std::endl;
         throw SGPP::base::operation_exception(errorString.str());
         }
         err = clSetKernelArg(kernel_multTrans[i], 2, sizeof(cl_mem),
         this->deviceMask.getBuffer(i));
         if (err != CL_SUCCESS) {
         std::stringstream errorString;
         errorString
         << "OCL Error: Failed to create kernel arguments for device "
         << i << std::endl;
         throw SGPP::base::operation_exception(errorString.str());
         }
         err = clSetKernelArg(kernel_multTrans[i], 3, sizeof(cl_mem),
         this->deviceOffset.getBuffer(i));
         if (err != CL_SUCCESS) {
         std::stringstream errorString;
         errorString
         << "OCL Error: Failed to create kernel arguments for device "
         << i << std::endl;
         throw SGPP::base::operation_exception(errorString.str());
         }
         err = clSetKernelArg(kernel_multTrans[i], 4, sizeof(cl_mem),
         this->deviceData.getBuffer(i));
         if (err != CL_SUCCESS) {
         std::stringstream errorString;
         errorString
         << "OCL Error: Failed to create kernel arguments for device "
         << i << std::endl;
         throw SGPP::base::operation_exception(errorString.str());
         }
         err = clSetKernelArg(kernel_multTrans[i], 5, sizeof(cl_mem),
         this->deviceTemp.getBuffer(i));
         if (err != CL_SUCCESS) {
         std::stringstream errorString;
         errorString
         << "OCL Error: Failed to create kernel arguments for device "
         << i << std::endl;
         throw SGPP::base::operation_exception(errorString.str());
         }
         err = clSetKernelArg(kernel_multTrans[i], 6, sizeof(cl_mem),
         this->deviceGrid.getBuffer(i));
         if (err != CL_SUCCESS) {
         std::stringstream errorString;
         errorString
         << "OCL Error: Failed to create kernel arguments for device "
         << i << std::endl;
         throw SGPP::base::operation_exception(errorString.str());
         }
         err = clSetKernelArg(kernel_multTrans[i], 7, sizeof(cl_uint),
         &clSourceSize); // sourceSize == number of entries in dataset
         if (err != CL_SUCCESS) {
         std::stringstream errorString;
         errorString
         << "OCL Error: Failed to create kernel arguments for device "
         << i << std::endl;
         throw SGPP::base::operation_exception(errorString.str());
         }
         err = clSetKernelArg(kernel_multTrans[i], 8, sizeof(cl_uint),
         &gpu_start_data);
         if (err != CL_SUCCESS) {
         std::stringstream errorString;
         errorString
         << "OCL Error: Failed to create kernel arguments for device "
         << i << std::endl;
         throw SGPP::base::operation_exception(errorString.str());
         }
         err = clSetKernelArg(kernel_multTrans[i], 9, sizeof(cl_uint),
         &gpu_end_data);
         if (err != CL_SUCCESS) {
         std::stringstream errorString;
         errorString
         << "OCL Error: Failed to create kernel arguments for device "
         << i << std::endl;
         throw SGPP::base::operation_exception(errorString.str());
         }
         }
         }

         cl_event* clTimings = new cl_event[manager->num_devices];

         // enqueue kernels
         size_t local = (*parameters)["LOCAL_SIZE"].getUInt();
         size_t active_devices = 0;

         for (size_t i = 0; i < manager->num_devices; i++) {
         size_t rangeSize = gpu_end_index_grid[i] - gpu_start_index_grid[i];

         if (rangeSize > 0) {
         //      std::cout << "enqueuing device: " << i << std::endl;
         err = clEnqueueNDRangeKernel(manager->command_queue[i], kernel_multTrans[i],
         1, &gpu_start_index_grid[i], &rangeSize, &local, 0, nullptr,
         &(clTimings[i]));

         if (active_devices != i) {
         std::stringstream errorString;
         errorString
         << "OCL Error: Multiple GPUs is erroneous, because only the last chunks may handle 0 entries. Synchronization will be not correct."
         << std::endl;
         throw SGPP::base::operation_exception(errorString.str());
         }

         active_devices++;

         if (err != CL_SUCCESS) {
         std::stringstream errorString;
         errorString
         << "OCL Error: Failed to enqueue kernel command! Error code: "
         << err << std::endl;
         throw SGPP::base::operation_exception(errorString.str());
         }
         }
         }

         deviceGrid.readFromBuffer(gpu_start_index_grid, gpu_end_index_grid);

         for (size_t i = start_index_grid; i < end_index_grid; i++) {
         result[i] = hostGrid[i];
         }

         // determine kernel execution time
         for (size_t i = 0; i < manager->num_devices; i++) {
         double tmpTime;
         cl_ulong startTime, endTime;
         startTime = endTime = 0;

         if (gpu_end_index_grid[i] > gpu_start_index_grid[i]) {
         err = clGetEventProfilingInfo(clTimings[i],
         CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime, nullptr);

         if (err != CL_SUCCESS) {
         std::stringstream errorString;
         errorString
         << "OCL Error: Failed to read start-time from command queue (or crash in multTranspose)! Error code: "
         << err << std::endl;
         throw SGPP::base::operation_exception(errorString.str());
         }

         err = clGetEventProfilingInfo(clTimings[i],
         CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, nullptr);

         if (err != CL_SUCCESS) {
         std::stringstream errorString;
         errorString
         << "OCL Error: Failed to read end-time from command queue! Error code: "
         << err << std::endl;
         throw SGPP::base::operation_exception(errorString.str());
         }
         }

         tmpTime = (double) (endTime - startTime);
         tmpTime *= 1e-9;
         this->deviceTimingsMultTranspose[i] = tmpTime;

         if (tmpTime > time) {
         time = tmpTime;
         }
         }

         // clean up
         for (size_t i = 0; i < manager->num_devices; i++) {
         if (gpu_end_index_grid[i] > gpu_start_index_grid[i]) {
         clReleaseEvent(clTimings[i]);
         }
         }

         delete[] gpu_start_index_grid;
         delete[] gpu_end_index_grid;

         delete[] clTimings;
         return time;
         */
        return 0.0;
    }
private:

    void releaseGridBuffers() {
        this->deviceLevel.freeBuffer();
        this->deviceIndex.freeBuffer();
        this->deviceMask.freeBuffer();
        this->deviceOffset.freeBuffer();
        this->deviceAlpha.freeBuffer();
    }

    void releaseDataBuffers() {
        this->deviceData.freeBuffer();
        this->deviceResultData.freeBuffer();
        this->deviceResultGrid.freeBuffer();
    }

    void initGridBuffers(std::vector<real_type> &level, std::vector<real_type> &index, std::vector<real_type> &mask,
            std::vector<real_type> &offset, std::vector<real_type> &alpha, size_t kernelStartGrid,
            size_t kernelEndGrid) {

        size_t gridRange = kernelEndGrid - kernelStartGrid;

        if (!this->deviceLevel.isInitialized()) {
            this->deviceLevel.initializeBuffer(gridRange * this->dims);
        }

        if (!this->deviceIndex.isInitialized()) {
            this->deviceIndex.initializeBuffer(gridRange * this->dims);
        }

        if (!this->deviceMask.isInitialized()) {
            this->deviceMask.initializeBuffer(gridRange * this->dims);
        }

        if (!this->deviceOffset.isInitialized()) {
            this->deviceOffset.initializeBuffer(gridRange * this->dims);
        }

        if (!this->deviceAlpha.isInitialized()) {
            this->deviceAlpha.initializeBuffer(gridRange);
        }

        std::vector<real_type> &deviceLevelHost = deviceLevel.getHostPointer();
        std::vector<real_type> &deviceIndexHost = deviceIndex.getHostPointer();
        std::vector<real_type> &deviceMaskHost = deviceMask.getHostPointer();
        std::vector<real_type> &deviceOffsetHost = deviceOffset.getHostPointer();
        std::vector<real_type> &deviceAlphaHost = deviceAlpha.getHostPointer();

        size_t deviceDataIndex = 0;
        for (size_t i = kernelStartGrid; i < kernelEndGrid; i++) {
            for (size_t d = 0; d < dims; d++) {
                deviceLevelHost[deviceDataIndex * dims + d] = level[i * dims + d];
                deviceIndexHost[deviceDataIndex * dims + d] = index[i * dims + d];
                deviceMaskHost[deviceDataIndex * dims + d] = mask[i * dims + d];
                deviceOffsetHost[deviceDataIndex * dims + d] = offset[i * dims + d];
            }
            deviceAlphaHost[deviceDataIndex] = alpha[i];
            deviceDataIndex += 1;
        }
        deviceLevel.writeToBuffer();
        deviceIndex.writeToBuffer();
        deviceMask.writeToBuffer();
        deviceOffset.writeToBuffer();
        deviceAlpha.writeToBuffer();
    }

    void initDatasetBuffers(std::vector<real_type> &dataset, size_t kernelStartData, size_t kernelEndData) {

        size_t range = kernelEndData - kernelStartData;

        if (!deviceData.isInitialized()) {
            this->deviceData.initializeBuffer(range * this->dims);
        } else if (range * this->dims != deviceData.size()) {
            deviceData.freeBuffer();
            this->deviceData.initializeBuffer(range * this->dims);
        }

        std::vector<real_type> &deviceDataHost = deviceData.getHostPointer();
        size_t dataPoints = dataset.size() / dims;
        for (size_t d = 0; d < dims; d++) {
            size_t deviceDataIndex = 0;
            for (size_t i = kernelStartData; i < kernelEndData; i++) {
                deviceDataHost[d * range + deviceDataIndex] = dataset[d * dataPoints + i];
                deviceDataIndex += 1;
            }
        }
        deviceData.writeToBuffer();
    }

    void initDatasetResultBuffers(size_t kernelStartData, size_t kernelEndData) {

        size_t range = kernelEndData - kernelStartData;

        if (!deviceResultData.isInitialized()) {
            deviceResultData.initializeBuffer(range);
        } else if (range != deviceResultData.size()) {
            deviceResultData.freeBuffer();
            deviceResultData.initializeBuffer(range);
            //all grid buffers can be reused, if already initialized
        }

        std::vector<real_type> &deviceResultDataHost = deviceResultData.getHostPointer();
        for (size_t i = 0; i < range; i++) {
            deviceResultDataHost[i] = 0.0;
        }
        deviceResultData.writeToBuffer();
    }

}
;

}
}
