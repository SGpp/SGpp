// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/opencl/OpenCLConfigurationParameters.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>

namespace SGPP {
namespace base {

OpenCLConfigurationParameters::OpenCLConfigurationParameters(
  std::string fileName,
  std::map<std::string, std::string> defaultParameters) {
  this->parameters["LOCAL_SIZE"] = "64";
  this->parameters["ENABLE_OPTIMIZATIONS"] = "true";
  this->parameters["OPTIMIZATION_FLAGS"] =
    "-cl-finite-math-only -cl-fast-relaxed-math";
  this->parameters["INTERNAL_PRECISION"] = "double";
  this->parameters["PLATFORM"] = "first";
  this->parameters["DEVICE_TYPE"] = "CL_DEVICE_TYPE_ALL";
  this->parameters["MAX_DEVICES"] = "0";
  this->parameters["SELECT_SPECIFIC_DEVICE"] = "DISABLED";
  this->parameters["REUSE_SOURCE"] = "false";
  this->parameters["OCL_MANAGER_VERBOSE"] = "false";

  this->readFromMap(defaultParameters);
  this->readFromFile(fileName);
}

}  // namespace base
}  // namespace SGPP
