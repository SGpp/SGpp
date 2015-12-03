#include <sgpp/base/opencl/OCLPlatformWrapper.hpp>

#include <sstream>

#include <sgpp/globaldef.hpp>

#include <sgpp/base/exception/operation_exception.hpp>

namespace SGPP {
namespace base {

OCLPlatformWrapper::OCLPlatformWrapper(cl_platform_id platformId, char (&platformName)[128], std::vector<cl_device_id> &deviceIds, std::vector<std::string> &deviceNames) :
        platformId(platformId), deviceIds(deviceIds), deviceNames(deviceNames) {

    for (size_t i = 0; i < 128; i++) {
        this->platformName[i] = platformName[i];
    }

    cl_int err;
    // Create OpenCL context
    cl_context_properties properties[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties) platformId, 0 };

    this->context = clCreateContext(properties, (cl_uint) this->deviceIds.size(),
            this->deviceIds.data(), nullptr, nullptr, &err);

    if (err != CL_SUCCESS) {
        std::stringstream errorString;
        errorString << "OCL Error: Failed to create OpenCL context! Error Code: " << err << std::endl;
        throw SGPP::base::operation_exception(errorString.str());
    }

    // Create a command queue for each device
    for (uint32_t i = 0; i < this->deviceIds.size(); i++) {
        this->commandQueues.push_back(
                clCreateCommandQueue(this->context, this->deviceIds[i],
                CL_QUEUE_PROFILING_ENABLE, &err));

        if (err != CL_SUCCESS) {
            std::stringstream errorString;
            errorString << "OCL Error: Failed to create command queue! Error Code: " << err << std::endl;
            throw SGPP::base::operation_exception(errorString.str());
        }
    }
}

size_t OCLPlatformWrapper::getDeviceCount() {
    return this->deviceIds.size();
}

}
}
