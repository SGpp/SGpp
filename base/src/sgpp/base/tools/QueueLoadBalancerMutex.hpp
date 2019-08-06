// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/globaldef.hpp>

#include <mutex>

namespace sgpp {
namespace base {

class QueueLoadBalancerMutex {
 private:
  bool isInitialized;
  size_t start;
  size_t end;
  size_t blockSize;
  size_t range;
  size_t currentStart;
  std::recursive_mutex nextSegmentMutex;

 public:
  // end is assumed to be padded for blocksize! might return segements that end
  // after end
  QueueLoadBalancerMutex()
      : isInitialized(false), start(0), end(0), blockSize(0), range(0), currentStart(0) {}

  void initialize(const size_t start, const size_t end, const size_t blockSize) {
    this->start = start;
    this->end = end;
    this->blockSize = blockSize;
    this->range = end - start;
    this->currentStart = start;
    this->isInitialized = true;
  }

  // is thread-safe
  bool getNextSegment(const size_t scheduleSize, size_t &segmentStart, size_t &segmentEnd) {
    if (!this->isInitialized) {
      throw base::operation_exception("QueueLoadBalancer: queue load balancer not initialized!");
    } else if (blockSize == 0) {
      throw base::operation_exception("QueueLoadBalancer: block size must not be zero!");
    } else if (scheduleSize % blockSize != 0) {
      throw base::operation_exception(
          "QueueLoadBalancer: schedule size is not divisible by block size!");
    } else if (start % blockSize != 0) {
      throw base::operation_exception("QueueLoadBalancer: start not divisible by block size");
    } else if (end % blockSize != 0) {
      throw base::operation_exception("QueueLoadBalancer: end not divisible by block size");
    }

    std::lock_guard<std::recursive_mutex> lock(nextSegmentMutex);

    bool segmentAvailable = true;
    if (currentStart == end) {
      segmentAvailable = false;
    } else {
      segmentStart = currentStart;
      if (currentStart + scheduleSize <= end) {
        segmentEnd = currentStart + scheduleSize;
      } else {
        segmentEnd = end;
      }

      currentStart = segmentEnd;
    }
    return segmentAvailable;
  }

  // is thread-safe
  void reset() {
    std::lock_guard<std::recursive_mutex> lock(nextSegmentMutex);
    currentStart = start;
  }
};
}  // namespace base
}  // namespace sgpp
