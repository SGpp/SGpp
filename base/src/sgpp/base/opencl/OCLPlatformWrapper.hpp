#pragma once

#include <vector>

#include <CL/cl.h>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace base {

class OCLPlatformWrapper {
 public:
  cl_platform_id platformId;
  char platformName[128];
  cl_context context;
  std::vector<cl_device_id> deviceIds;
  std::vector<cl_command_queue> commandQueues;

  OCLPlatformWrapper(cl_platform_id platformId, char (&platformName)[128],
                     std::vector<cl_device_id>& deviceIds);

  size_t getDeviceCount();
};

}  // namespace base
}  // namespace SGPP
