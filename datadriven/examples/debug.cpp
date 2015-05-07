/*
 * multiEvalPerformance.cpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include <iostream>
#include <random>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLStretchedBufferMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLClonedBufferMultiPlatform.hpp>

int main(int argc, char** argv) {

    std::map<std::string, std::string> defaultParameter;
    defaultParameter["KERNEL_USE_LOCAL_MEMORY"] = "true";
    defaultParameter["KERNEL_DATA_BLOCKING_SIZE"] = "1";
    defaultParameter["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";
    //  defaultParameter["KERNEL_GRID_BLOCK_SIZE"] = "1";
    defaultParameter["KERNEL_TRANS_DATA_BLOCK_SIZE"] = "1";
    defaultParameter["KERNEL_TRANS_UNROLL_1D"] = "true";
    defaultParameter["KERNEL_STORE_DATA"] = "array";

    SGPP::base::OpenCLConfigurationParameters parameters("StreamingModOCLFastMultiPlatform.cfg", defaultParameter);
    SGPP::base::OCLManagerMultiPlatform manager(parameters);

    SGPP::base::OCLClonedBufferMultiPlatform clonedBuffer(manager);

    double *clonedBufferHost = new double[100];
    clonedBuffer.initializeBuffer(clonedBufferHost, sizeof(double), 100);

    clonedBuffer.writeToBuffer(clonedBufferHost);

    clonedBuffer.readFromBuffer(clonedBufferHost);


    SGPP::base::OCLStretchedBufferMultiPlatform buffer(manager);

    buffer.initializeBuffer(sizeof(double), 100);

    double *hostBuffer = (double *) buffer.getMappedHostBuffer(manager.platforms[0].platformId);

    for (size_t i = 0; i < 100; i++) {
        hostBuffer[i] = static_cast<double>(i);
    }

    std::cout << "copying to other buffers" << std::endl;
    buffer.copyToOtherHostBuffers(manager.platforms[0].platformId);
    std::cout << "write to buffer" << std::endl;
    buffer.writeToBuffer();

    std::map<cl_platform_id, size_t *> indexStart;
    for (SGPP::base::OCLPlatformWrapper &platform : manager.platforms) {
        indexStart[platform.platformId] = new size_t[platform.deviceCount];
        for (size_t i = 0; i < platform.deviceCount; i++) {
            indexStart[platform.platformId][i] = i;
        }
    }

    std::map<cl_platform_id, size_t *> indexEnd;
    for (SGPP::base::OCLPlatformWrapper &platform : manager.platforms) {
        indexEnd[platform.platformId] = new size_t[platform.deviceCount];
        for (size_t i = 0; i < platform.deviceCount; i++) {
            indexEnd[platform.platformId][i] = i + 1;
        }
    }

    std::cout << "read from buffer" << std::endl;
    buffer.readFromBuffer(indexStart, indexEnd);
    std::cout << "combine buffers in single host buffer" << std::endl;
    buffer.combineBuffer(indexStart, indexEnd, manager.platforms[0].platformId);
    for (size_t i = 0; i < 100; i++) {
        if (i != 0) {
            std::cout  << ", ";
        }
        std::cout << hostBuffer[i];
    }
    std::cout << std::endl;

    std::cout << "all done" << std::endl;

}

