/*
 * LinearLoadBalancer.hpp
 *
 *  Created on: Mar 30, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace base {

class QueueLoadBalancer {
private:
    bool isInitialized;
    size_t start;
    size_t end;
    size_t range;
    size_t currentStart;

public:
    //end is assumed to be padded for blocksize! might return segements that end after end
    QueueLoadBalancer() :
            isInitialized(false), start(0), end(0), range(0), currentStart(0) {
    }

    void initialize(const size_t start, const size_t end) {
        this->start = start;
        this->end = end;
        this->range = end - start;
        this->currentStart = start;
        this->isInitialized = true;
    }

    //is thread-safe
    bool getNextSegment(const size_t scheduleSize, const size_t blockSize, size_t &segmentStart, size_t &segmentEnd) {

        if (!this->isInitialized) {
            throw base::operation_exception("QueueLoadBalancer: queue load balancer not initialized!");
        } else if (blockSize == 0) {
            throw base::operation_exception("QueueLoadBalancer: block size must not be zero!");
        } else if (scheduleSize % blockSize != 0) {
            throw base::operation_exception("QueueLoadBalancer: schedule size is not divisible by block size!");
        } else if (start % blockSize != 0) {
            throw base::operation_exception("QueueLoadBalancer: start not divisible by block size");
        } else if (end % blockSize != 0) {
            throw base::operation_exception("QueueLoadBalancer: end not divisible by block size");
        }

        bool segmentAvailable = true;
#pragma omp critical (QueueLoadBalancer_getNextSegment)
        {
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
        }
        return segmentAvailable;
    }

    // is thread-safe
    void reset() {
#pragma omp critical (QueueLoadBalancer_getNextSegment)
        {
            currentStart = start;
        }
    }

}
;

}
}

