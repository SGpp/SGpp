/*
 * multiEvalPerformance.cpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include <iostream>
#if USE_OCL == 1

#include <random>
#include <thread>
#include <chrono>

#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLClonedBufferMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLStretchedBufferMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLConfigurationParameters.hpp>

#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/base/opencl/OCLClonedBuffer.hpp>
#include <sgpp/base/opencl/OCLStretchedBuffer.hpp>

using namespace SGPP::base;

void doStuffOld(std::shared_ptr<OCLManager> manager, double *values, size_t valueSize) {
//    OCLClonedBuffer buffer(manager);
//    buffer.initializeBuffer(values, sizeof(double), valueSize);

    OCLStretchedBuffer stretched(manager);
    stretched.initializeBuffer(sizeof(double), valueSize);
}

void doStuff(std::shared_ptr<OCLManagerMultiPlatform> manager, double *values, size_t valueSize) {
//    OCLClonedBufferMultiPlatform buffer(manager);
//    buffer.initializeBuffer(values, sizeof(double), valueSize);

    OCLStretchedBufferMultiPlatform stretched(manager);
    stretched.freeBuffer();
    stretched.initializeBuffer(sizeof(double), valueSize);
    stretched.freeBuffer();
    stretched.initializeBuffer(sizeof(double), valueSize);
}

int main(int argc, char** argv) {

    auto parameters = std::make_shared<OCLConfigurationParameters>();

    auto manager = std::make_shared<OCLManagerMultiPlatform>(parameters);

    const size_t valueSize = 100000000;

    double *values = new double[valueSize];
    for (size_t i = 0; i < valueSize; i++) {
        values[i] = static_cast<double>(i);
    }

    const size_t iterations = 10000;

    auto managerOld = std::make_shared<OCLManager>(parameters);

    for (size_t i = 0; i < iterations; i++) {
        std::cout << "it: " << i << std::endl;
        doStuff(manager, values, valueSize);
//        doStuffOld(managerOld, values, valueSize);
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

}
#else
#include <iostream>
int main(int argc, char** argv) {
    std::cout << "This examples requires OpenCL to be enabled. (build with USE_OCL=1)" << std::endl;
	return 0;
}
#endif
