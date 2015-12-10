#pragma once

#include <vector>

#include <CL/cl.h>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace base {

class OCLDevice {
public:
    cl_platform_id platformId;
    cl_device_id deviceId;
    std::string platformName;
    std::string deviceName;
    cl_context context;
    cl_command_queue commandQueue;

    OCLDevice(cl_platform_id platformId, cl_device_id deviceId, const std::string &platformName,
            const std::string &deviceName, cl_context context, cl_command_queue commandQueue) :
            platformId(platformId), deviceId(deviceId), platformName(platformName), deviceName(deviceName), context(
                    context), commandQueue(commandQueue) {

    }
};

}
}
