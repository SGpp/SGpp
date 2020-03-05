// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/opencl/OCLClonedBufferMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLStretchedBufferMultiPlatform.hpp>

#include <sgpp/base/opencl/OCLClonedBuffer.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/base/opencl/OCLStretchedBuffer.hpp>

#include <iostream>

#include <chrono>
#include <random>
#include <thread>  // NOLINT(build/c++11)


void doStuffOld(std::shared_ptr<sgpp::base::OCLManager> manager, double* values, size_t valueSize) {
  //    OCLClonedBuffer buffer(manager);
  //    buffer.initializeBuffer(values, sizeof(double), valueSize);

  sgpp::base::OCLStretchedBuffer stretched(manager);
  stretched.initializeBuffer(sizeof(double), valueSize);
}

void doStuff(std::shared_ptr<sgpp::base::OCLManagerMultiPlatform> manager, double* values,
             size_t valueSize) {
  //    OCLClonedBufferMultiPlatform buffer(manager);
  //    buffer.initializeBuffer(values, sizeof(double), valueSize);

  sgpp::base::OCLStretchedBufferMultiPlatform stretched(manager);
  stretched.freeBuffer();
  stretched.initializeBuffer(sizeof(double), valueSize);
  stretched.freeBuffer();
  stretched.initializeBuffer(sizeof(double), valueSize);
}

int main(int argc, char** argv) {
  auto parameters = std::make_shared<sgpp::base::OCLOperationConfiguration>();

  auto manager = std::make_shared<sgpp::base::OCLManagerMultiPlatform>(parameters);

  const size_t valueSize = 100000000;

  double* values = new double[valueSize];

  for (size_t i = 0; i < valueSize; i++) {
    values[i] = static_cast<double>(i);
  }

  const size_t iterations = 10000;

  auto managerOld = std::make_shared<sgpp::base::OCLManager>(parameters);

  for (size_t i = 0; i < iterations; i++) {
    std::cout << "it: " << i << std::endl;
    doStuff(manager, values, valueSize);
    //        doStuffOld(managerOld, values, valueSize);
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }
}
