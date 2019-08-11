// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <vector>
#include <map>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/tools/ConfigurationParameters.hpp>
#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>

namespace sgpp {
namespace base {

class LinearLoadBalancerMultiPlatform {
 private:
  std::shared_ptr<OCLManagerMultiPlatform> manager;
  std::shared_ptr<base::OCLOperationConfiguration> parameters;
  std::map<cl_platform_id, std::vector<double>> weights;
  std::map<cl_platform_id, std::vector<double>> partition;
  std::map<cl_platform_id, std::vector<double>> lastMeaningfulPartition;

 public:
  LinearLoadBalancerMultiPlatform(std::shared_ptr<OCLManagerMultiPlatform> manager,
                                  std::shared_ptr<base::OCLOperationConfiguration> parameters)
      : manager(manager), parameters(parameters) {
    for (OCLPlatformWrapper &platform : manager->platforms) {
      this->weights[platform.platformId] = std::vector<double>(platform.getDeviceCount());
      this->partition[platform.platformId] = std::vector<double>(platform.getDeviceCount());
      this->lastMeaningfulPartition[platform.platformId] =
          std::vector<double>(platform.getDeviceCount());
    }

    for (OCLPlatformWrapper &platform : manager->platforms) {
      for (size_t i = 0; i < platform.getDeviceCount(); i++) {
        // initialize with same timing to
        // enforce equal problem sizes in the beginning
        this->partition[platform.platformId][i] =
            1.0 / static_cast<double>(manager->overallDeviceCount);
        this->lastMeaningfulPartition[platform.platformId][i] =
            1.0 / static_cast<double>(manager->overallDeviceCount);
      }
    }
  }

  void getPartitionSegments(size_t start, size_t end, size_t blockSize,
                            std::map<cl_platform_id, size_t *> segmentStart,
                            std::map<cl_platform_id, size_t *> segmentEnd) {
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

    for (OCLPlatformWrapper &platform : manager->platforms) {
      for (size_t i = 0; i < platform.getDeviceCount(); i++) {
        size_t partitionElements =
            static_cast<size_t>(static_cast<double>(totalSize) * partition[platform.platformId][i]);

        if (currentStartIndex != end && partitionElements == 0) {
          partitionElements = 1;
        }

        // last device has to ensure that all data is in one partition
        if (currentStartIndex + partitionElements > end ||
            ((&platform == &(manager->platforms[manager->platforms.size() - 1])) &&
             (i == platform.getDeviceCount() - 1))) {
          partitionElements = end - currentStartIndex;
        }

        // employ padding
        size_t remainder = partitionElements % blockSize;
        size_t padding = 0;

        if (remainder != 0) {
          padding = blockSize - remainder;
        }

        partitionElements += padding;

        segmentStart[platform.platformId][i] = currentStartIndex;
        segmentEnd[platform.platformId][i] = currentStartIndex + partitionElements;

        if (setVerboseLoadBalancing) {
          std::cout << "device: " << i << " from: " << segmentStart[platform.platformId][i]
                    << " to: " << segmentEnd[platform.platformId][i] << std::endl;
        }

        currentStartIndex += partitionElements;
      }
    }
  }

  void update(std::map<cl_platform_id, double *> timings) {
    bool setVerboseLoadBalancing = (*parameters)["LOAD_BALANCING_VERBOSE"].getBool();

    // recalculate weights
    for (OCLPlatformWrapper platform : manager->platforms) {
      for (size_t i = 0; i < platform.getDeviceCount(); i++) {
        weights[platform.platformId][i] =
            timings[platform.platformId][i] / partition[platform.platformId][i];

        if (setVerboseLoadBalancing) {
          std::cout << "platform: \"" << platform.platformName << "\" device: " << i << " took "
                    << timings[platform.platformId][i] << "s" << std::endl;
        }
      }
    }

    // calculate the optimal duration
    double t = 0.0;

    for (OCLPlatformWrapper platform : manager->platforms) {
      for (size_t i = 0; i < platform.getDeviceCount(); i++) {
        t += 1.0 / weights[platform.platformId][i];
      }
    }

    t = 1.0 / t;

    if (setVerboseLoadBalancing) {
      std::cout << "t: " << t << std::endl;
    }

    // calculate optimal partition
    for (OCLPlatformWrapper platform : manager->platforms) {
      for (size_t i = 0; i < platform.getDeviceCount(); i++) {
        if (t == 0.0) {
          // partition[platform.platformId][i] =
          // 1.0 / static_cast<double>(manager->overallDeviceCount);
          partition[platform.platformId][i] = lastMeaningfulPartition[platform.platformId][i];
        } else {
          partition[platform.platformId][i] = t / weights[platform.platformId][i];
          lastMeaningfulPartition[platform.platformId][i] = partition[platform.platformId][i];
        }

        if (setVerboseLoadBalancing) {
          std::cout << "device: " << i << " partition size: " << partition[platform.platformId][i]
                    << std::endl;
        }
      }
    }
  }
};

}  // namespace base
}  // namespace sgpp
