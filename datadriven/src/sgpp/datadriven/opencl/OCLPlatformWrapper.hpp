#pragma once

#include <CL/cl.h>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace base {

class OCLPlatformWrapper {
public:
    cl_platform_id platformId;
    char platformName[128];
    cl_context context;
    cl_device_id *deviceIds;
    size_t deviceCount;
    cl_command_queue *commandQueues;

    OCLPlatformWrapper(cl_platform_id platformId);

    OCLPlatformWrapper(const OCLPlatformWrapper &original);

    ~OCLPlatformWrapper();
};

}
}
