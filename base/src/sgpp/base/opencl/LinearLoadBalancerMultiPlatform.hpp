/*
 * LinearLoadBalancer.hpp
 *
 *  Created on: Mar 30, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/ConfigurationParameters.hpp>

#include "OCLManagerMultiPlatform.hpp"

namespace SGPP {
namespace base {

class LinearLoadBalancerMultiPlatform {
private:
    OCLManagerMultiPlatform &manager;
    base::OCLConfigurationParameters& parameters;
    std::map<cl_platform_id, double *> weights;
    std::map<cl_platform_id, double *> partition;
public:
    LinearLoadBalancerMultiPlatform(OCLManagerMultiPlatform& manager, base::OCLConfigurationParameters& parameters) :
            manager(manager), parameters(parameters) {
        for (OCLPlatformWrapper platform : manager.platforms) {
            this->weights[platform.platformId] = new double[platform.deviceCount];
            this->partition[platform.platformId] = new double[platform.deviceCount];
        }

        for (OCLPlatformWrapper platform : manager.platforms) {
            for (size_t i = 0; i < platform.deviceCount; i++) {
                //initialize with same timing to enforce equal problem sizes in the beginning
                this->partition[platform.platformId][i] = 1.0 / static_cast<double>(manager.overallDeviceCount);
            }
        }
    }

    ~LinearLoadBalancerMultiPlatform() {
        for (OCLPlatformWrapper platform : manager.platforms) {
            delete[] this->weights[platform.platformId];
            delete[] this->partition[platform.platformId];
        }
    }

    void getPartitionSegments(size_t start, size_t end, size_t blockSize,
            std::map<cl_platform_id, size_t *> segmentStart, std::map<cl_platform_id, size_t *> segmentEnd) {
        bool setVerboseLoadBalancing = parameters.getAsBoolean("LOAD_BALANCING_VERBOSE");
        size_t totalSize = end - start;

        // check for valid input
        if (blockSize == 0) {
            throw SGPP::base::operation_exception("blockSize must not be zero!");
        }

        if (totalSize % blockSize != 0) {
            throw SGPP::base::operation_exception(
                    "totalSize must be divisible by blockSize without remainder, but it is not!");
        }

        size_t currentStartIndex = start;

        for (OCLPlatformWrapper &platform : manager.platforms) {
            for (size_t i = 0; i < platform.deviceCount; i++) {
                size_t partitionElements = static_cast<size_t>(static_cast<double>(totalSize)
                        * partition[platform.platformId][i]);

                if (currentStartIndex != end && partitionElements == 0) {
                    partitionElements = 1;
                }

                //last device has to ensure that all data is in one partition
                if (currentStartIndex + partitionElements > end
                        || ((&platform == &(manager.platforms[manager.platforms.size() - 1]))
                                && (i == platform.deviceCount - 1))) {
                    partitionElements = end - currentStartIndex;
                }

                //employ padding
                size_t remainder = partitionElements % blockSize;
                size_t padding = 0;

                if (remainder != 0) {
                    padding = blockSize - remainder;
                }

                partitionElements += padding;

                segmentStart[platform.platformId][i] = currentStartIndex;
                segmentEnd[platform.platformId][i] = currentStartIndex + partitionElements;

                if (setVerboseLoadBalancing) {
                    std::cout << "device: " << i << " from: " << segmentStart[platform.platformId][i] << " to: "
                            << segmentEnd[platform.platformId][i] << std::endl;
                }

                currentStartIndex += partitionElements;
            }
        }
    }

    //TODO: consider inactive device due to nothing to do?
    void update(std::map<cl_platform_id, double *> timings) {
        bool setVerboseLoadBalancing = parameters.getAsBoolean("LOAD_BALANCING_VERBOSE");

        //recalculate weights
        for (OCLPlatformWrapper platform : manager.platforms) {
            for (size_t i = 0; i < platform.deviceCount; i++) {
                weights[platform.platformId][i] = timings[platform.platformId][i] / partition[platform.platformId][i];
            }
        }

        //calculate the optimal duration
        double t = 0.0;

        for (OCLPlatformWrapper platform : manager.platforms) {
            for (size_t i = 0; i < platform.deviceCount; i++) {
                t += 1.0 / weights[platform.platformId][i];
            }
        }

        t = 1.0 / t;

        if (setVerboseLoadBalancing) {
            std::cout << "t: " << t << std::endl;
        }

        //calculate optimal partition
        for (OCLPlatformWrapper platform : manager.platforms) {
            for (size_t i = 0; i < platform.deviceCount; i++) {
                if (t == 0.0) {
                    partition[platform.platformId][i] = 1.0 / static_cast<double>(manager.overallDeviceCount);
                } else {
                    partition[platform.platformId][i] = t / weights[platform.platformId][i];
                }

                if (setVerboseLoadBalancing) {
                    std::cout << "device: " << i << " partition size: " << partition[platform.platformId][i]
                            << std::endl;
                }
            }
        }
    }

}
;

}
}

