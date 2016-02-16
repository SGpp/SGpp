// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <vector>
#include <string>

#include "CL/cl.h"
#include "sgpp/globaldef.hpp"

namespace SGPP {
namespace base {

class OCLPlatformWrapper {
 public:
  cl_platform_id platformId;
  char platformName[128];
  cl_context context;
  std::vector<cl_device_id> deviceIds;
  std::vector<std::string> deviceNames;
  std::vector<cl_command_queue> commandQueues;

  OCLPlatformWrapper(cl_platform_id platformId, char(&platformName)[128],
                     const std::vector<cl_device_id> &deviceIds,
                     const std::vector<std::string> &deviceName);

  size_t getDeviceCount();
};
}  // namespace base
}  // namespace SGPP
