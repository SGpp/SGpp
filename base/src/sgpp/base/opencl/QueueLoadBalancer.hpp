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
    const size_t scheduleSize; //excl. blocksize
    size_t paddedScheduleSize;
    const size_t start;
    const size_t end;
    const size_t blockSize;
    const size_t range;
    size_t currentStart;

public:
    //end is assumed to be padded for blocksize! might return segements that end after end
    QueueLoadBalancer(const size_t scheduleSize, const size_t start, const size_t end, const size_t blockSize) :
            scheduleSize(scheduleSize), start(start), end(end), blockSize(blockSize), range(end - start), currentStart(
                    start) {

        // check for valid input
        if (blockSize == 0) {
            throw SGPP::base::operation_exception("blockSize must not be zero!");
        }

//        if (range % blockSize != 0) {
//            throw SGPP::base::operation_exception(
//                    "totalSize must be divisible by blockSize without remainder, but it is not!");
//        }

        //employ padding
        size_t remainder = scheduleSize % blockSize;

        if (remainder != 0) {
            paddedScheduleSize = scheduleSize + (blockSize - remainder);
        }
    }

    //is thread-safe
    bool getNextSegment(size_t &segmentStart, size_t &segmentEnd) {
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

