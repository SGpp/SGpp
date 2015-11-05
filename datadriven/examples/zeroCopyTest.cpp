/*
 * zeroCopyTest.cpp
 *
 *  Created on: Oct 21, 2015
 *      Author: leiterrl
 */

#include <iostream>
#if USE_OCL == 1
#include <random>
#include <thread>
#include <chrono>

#include <sgpp/base/opencl/OCLConfigurationParameters.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/base/opencl/OCLClonedBuffer.hpp>
#include <sgpp/base/opencl/OCLStretchedBuffer.hpp>
#include <sgpp/base/opencl/OCLZeroCopyBuffer.hpp>

using namespace SGPP::base;
using namespace std::chrono;

void testCloned(std::shared_ptr<OCLManager> manager, double *values, size_t valueSize) {
    OCLClonedBuffer cloned(manager);
    cloned.initializeBuffer(values, sizeof(double), valueSize);
}

void testStretched(std::shared_ptr<OCLManager> manager, double *values, size_t valueSize) {
    OCLStretchedBuffer stretched(manager);
    stretched.initializeBuffer(sizeof(double), valueSize);
    double* mappedBuffer = (double*)stretched.getMappedHostBuffer();

    for (size_t i = 0; i < valueSize; i++)
    {
        mappedBuffer[i] = values[i];
    }
}

void testZeroCopyReadOnly(std::shared_ptr<OCLManager> manager, double *values, size_t valueSize) {
    OCLZeroCopyBuffer buffer(manager);
    buffer.initializeBuffer( values, sizeof(double), valueSize, true);
}

int main(int argc, char** argv) {

    auto parameters = std::make_shared<OCLConfigurationParameters>();
    //parameters->set("OCL_MANAGER_VERBOSE", "true");
    parameters->set("MAX_DEVICES", "1");
    parameters->set("SELECT_SPECIFIC_DEVICE", "0");

    auto manager = std::make_shared<OCLManager>(parameters);

    const size_t valueSize = 25264128; //max buffer size for Spectre

    double *values = new double[valueSize];
    for (size_t i = 0; i < valueSize; i++) {
        values[i] = static_cast<double>(i);
    }

    const size_t iterations = 100;

    size_t sum = 0;
    double avg = 0.0;

    std::cout << "Cloned Buffer" << std::endl;
    for (size_t i = 0; i < iterations; i++) {
        std::cout << ".";
        std::cout.flush();

        high_resolution_clock::time_point t1 = high_resolution_clock::now();

        testCloned(manager, values, valueSize);

        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
        //std::cout << "duration: " << duration << std::endl;
        sum += duration;

        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    avg = (static_cast<double>(sum)/static_cast<double>(iterations)) * 0.001 * 0.001;
    std::cout << "Average: " << avg << std::endl;
    sum = 0;
    std::cout << std::endl;

    std::cout << "Stretched Buffer" << std::endl;
    for (size_t i = 0; i < iterations; i++) {
        std::cout << ".";
        std::cout.flush();

        high_resolution_clock::time_point t1 = high_resolution_clock::now();

        testStretched(manager,values,valueSize);

        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
        //std::cout << "duration: " << duration << std::endl;
        sum += duration;

        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    avg = (static_cast<double>(sum)/static_cast<double>(iterations)) * 0.001 * 0.001;
    std::cout << "Average: " << avg << std::endl;
    sum = 0;
    std::cout << std::endl;

    std::cout << "ZeroCopy Buffer" << std::endl;
    for (size_t i = 0; i < iterations; i++) {
        std::cout << ".";
        std::cout.flush();

        high_resolution_clock::time_point t1 = high_resolution_clock::now();

        testZeroCopyReadOnly(manager,values,valueSize);

        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
        //std::cout << "duration: " << duration << std::endl;
        sum += duration;

        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    avg = (static_cast<double>(sum)/static_cast<double>(iterations)) * 0.001 * 0.001;
    std::cout << "Average: " << avg << std::endl;
    sum = 0;
}
#else
int main(int argc, char** argv) {
    std::cout << "This examples requires OpenCL to be enabled. (build with USE_OCL=1)" << std::endl;
        return 0;
}
#endif

