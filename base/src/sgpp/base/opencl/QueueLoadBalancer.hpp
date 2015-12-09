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
    size_t scheduleSize; //excl. blocksize
    size_t paddedScheduleSize;
    size_t start;
    size_t end;
    size_t blockSize;
    size_t range;
    size_t currentStart;

public:
    //end is assumed to be padded for blocksize! might return segements that end after end
    QueueLoadBalancer() :
            isInitialized(false), scheduleSize(0), paddedScheduleSize(0), start(0), end(0), blockSize(0), range(0), currentStart(0) {
    }

    void initialize(const size_t scheduleSize, const size_t start, const size_t end, const size_t blockSize) {
        // check for valid input
        if (blockSize == 0) {
            throw SGPP::base::operation_exception("QueueLoadBalancer: blockSize must not be zero!");
        }

        if (start % blockSize != 0) {
            throw SGPP::base::operation_exception("QueueLoadBalancer: start not divisible by blockSize");
        }

        if (end % blockSize != 0) {
            throw SGPP::base::operation_exception("QueueLoadBalancer: end not divisible by blockSize");
        }

        this->scheduleSize = scheduleSize;
        this->start = start;
        this->end = end;
        this->range = end - start;
        this->blockSize = blockSize;
        this->currentStart = start;

        //employ padding
        size_t remainder = this->scheduleSize % blockSize;

        if (remainder != 0) {
            paddedScheduleSize = this->scheduleSize + (blockSize - remainder);
        } else {
            paddedScheduleSize = this->scheduleSize;
        }
        this->isInitialized = true;
    }

    //is thread-safe
    bool getNextSegment(size_t &segmentStart, size_t &segmentEnd) {

        if (!this->isInitialized) {
            throw SGPP::base::operation_exception("QueueLoadBalancer: queue load balancer not initialized!");
        }

        bool segmentAvailable = true;
#pragma omp critical (QueueLoadBalancer_getNextSegment)
        {

            if (currentStart == end) {
                segmentAvailable = false;
            } else {

                segmentStart = currentStart;
                if (currentStart + paddedScheduleSize <= end) {
                    segmentEnd = currentStart + paddedScheduleSize;
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

