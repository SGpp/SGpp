// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <CL/cl.h>
#include <vector>
#include <string>

#include <sgpp/globaldef.hpp>

namespace sgpp {
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
            const std::string &deviceName, cl_context context, cl_command_queue commandQueue)
      : platformId(platformId),
        deviceId(deviceId),
        platformName(platformName),
        deviceName(deviceName),
        context(context),
        commandQueue(commandQueue) {}
};
}  // namespace base
}  // namespace sgpp
