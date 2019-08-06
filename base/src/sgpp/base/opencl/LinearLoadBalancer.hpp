// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

namespace sgpp {
namespace base {

class LinearLoadBalancer {
 private:
  size_t deviceCount;
  std::shared_ptr<base::OCLOperationConfiguration> parameters;
  double* weights;
  double* partition;

 public:
  LinearLoadBalancer(std::shared_ptr<OCLManager> manager,
                     std::shared_ptr<base::OCLOperationConfiguration> parameters)
      : deviceCount(manager->num_devices), parameters(parameters) {
    this->weights = new double[manager->num_devices];
    this->partition = new double[manager->num_devices];

    for (size_t i = 0; i < manager->num_devices; i++) {
      // initialize with same timing to
      // enforce equal problem sizes in the beginning
      this->partition[i] = 1.0 / static_cast<double>(manager->num_devices);
    }
  }

  ~LinearLoadBalancer() {
    delete[] this->weights;
    delete[] this->partition;
  }

  void getPartitionSegments(size_t start, size_t end, size_t blockSize, size_t* segmentStart,
                            size_t* segmentEnd) {
    bool setVerboseLoadBalancing = (*parameters)["LOAD_BALANCING_VERBOSE"].getBool();
    size_t totalSize = end - start;

    // check for valid input
    if (blockSize == 0) {
      throw sgpp::base::operation_exception("blockSize must not be zero!");
    }

    if (totalSize % blockSize != 0) {
      throw sgpp::base::operation_exception(
          "totalSize must be divisible by blockSize without remainder, "
          "but it is not!");
    }

    size_t currentStartIndex = start;

    for (size_t i = 0; i < this->deviceCount; i++) {
      size_t partitionElements = static_cast<size_t>(static_cast<double>(totalSize) * partition[i]);

      if (currentStartIndex != end && partitionElements == 0) {
        partitionElements = 1;
      }

      // last device has to ensure that all data is in one partition
      if (currentStartIndex + partitionElements > end || i == this->deviceCount - 1) {
        partitionElements = end - currentStartIndex;
      }

      // employ padding
      size_t remainder = partitionElements % blockSize;
      size_t padding = 0;

      if (remainder != 0) {
        padding = blockSize - remainder;
      }

      partitionElements += padding;

      segmentStart[i] = currentStartIndex;
      segmentEnd[i] = currentStartIndex + partitionElements;

      if (setVerboseLoadBalancing) {
        std::cout << "device: " << i << " from: " << segmentStart[i] << " to: " << segmentEnd[i]
                  << std::endl;
      }

      currentStartIndex += partitionElements;
    }
  }

  void update(double* timings) {
    bool setVerboseLoadBalancing = (*parameters)["LOAD_BALANCING_VERBOSE"].getBool();

    // recalculate weights
    for (size_t i = 0; i < this->deviceCount; i++) {
      weights[i] = timings[i] / partition[i];
    }

    // calculate the optimal duration
    double t = 0.0;

    for (size_t i = 0; i < this->deviceCount; i++) {
      t += 1.0 / weights[i];
    }

    t = 1.0 / t;

    if (setVerboseLoadBalancing) {
      std::cout << "t: " << t << std::endl;
    }

    // calculate optimal partition
    for (size_t i = 0; i < this->deviceCount; i++) {
      if (t == 0.0) {
        partition[i] = 1.0 / static_cast<double>(this->deviceCount);
      } else {
        partition[i] = t / weights[i];
      }

      if (setVerboseLoadBalancing) {
        std::cout << "device: " << i << " partition size: " << partition[i] << std::endl;
      }
    }
  }
};

}  // namespace base
}  // namespace sgpp
