#include "../../../../../datadriven/src/sgpp/datadriven/opencl/OCLPlatformWrapper.hpp"

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace base {

OCLPlatformWrapper::OCLPlatformWrapper(cl_platform_id platformId) :
        platformId(platformId), context(nullptr), deviceIds(nullptr), deviceCount(0), commandQueues(nullptr) {
    for (size_t i = 0; i < 128; i++) {
        platformName[i] = '\0';
    }
}

OCLPlatformWrapper::OCLPlatformWrapper(const OCLPlatformWrapper &original) {
    platformId = original.platformId;

    for (size_t i = 0; i < 128; i++) {
        platformName[i] = original.platformName[i];
    }

    context = original.context;
    deviceCount = original.deviceCount;
    if (deviceCount > 0) {
        if (original.deviceIds != nullptr) {
            deviceIds = new cl_device_id[deviceCount];
            for (size_t i = 0; i < deviceCount; i++) {
                deviceIds[i] = original.deviceIds[i];
            }
        } else {
            deviceIds = nullptr;
        }
        if (original.commandQueues != nullptr) {
            commandQueues = new cl_command_queue[deviceCount];
            for (size_t i = 0; i < deviceCount; i++) {
                commandQueues[i] = original.commandQueues[i];
            }
        } else {
            commandQueues = nullptr;
        }
    } else {
        deviceIds = nullptr;
        commandQueues = nullptr;
    }
}

OCLPlatformWrapper::~OCLPlatformWrapper() {
    if (deviceIds != nullptr) {
        delete[] deviceIds;
    }
    if (commandQueues != nullptr) {
        delete[] commandQueues;
    }
}

}
}
