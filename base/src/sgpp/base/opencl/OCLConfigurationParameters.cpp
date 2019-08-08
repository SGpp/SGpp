// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <string>

#include <sgpp/base/opencl/OCLConfigurationParameters.hpp>

namespace sgpp {
namespace base {

OCLConfigurationParameters::OCLConfigurationParameters(
    std::string fileName,
    std::map<std::string, std::string> defaultParameters) {
  this->parameters["LOCAL_SIZE"] = "64";
  this->parameters["ENABLE_OPTIMIZATIONS"] = "true";
  this->parameters["OPTIMIZATION_FLAGS"] = "";
  this->parameters["INTERNAL_PRECISION"] = "double";
  this->parameters["PLATFORM"] = "first";
  this->parameters["DEVICE_TYPE"] = "CL_DEVICE_TYPE_ALL";
  this->parameters["MAX_DEVICES"] = "DISABLED";
  this->parameters["SELECT_SPECIFIC_DEVICE"] = "DISABLED";
  this->parameters["REUSE_SOURCE"] = "false";
  this->parameters["WRITE_SOURCE"] = "false";
  // sets the kernel to verbose
  this->parameters["VERBOSE"] = "true";
  // sets the manager to verbose
  this->parameters["OCL_MANAGER_VERBOSE"] = "false";
  this->parameters["LOAD_BALANCING_VERBOSE"] = "false";
  this->parameters["SHOW_BUILD_LOG"] = "false";

  this->readFromMap(defaultParameters);
  this->readFromFile(fileName);
}

OCLConfigurationParameters::OCLConfigurationParameters() {
  this->parameters["LOCAL_SIZE"] = "64";
  this->parameters["ENABLE_OPTIMIZATIONS"] = "true";
  this->parameters["OPTIMIZATION_FLAGS"] = "";
  this->parameters["INTERNAL_PRECISION"] = "double";
  this->parameters["PLATFORM"] = "first";
  this->parameters["DEVICE_TYPE"] = "CL_DEVICE_TYPE_ALL";
  this->parameters["MAX_DEVICES"] = "DISABLED";
  this->parameters["SELECT_SPECIFIC_DEVICE"] = "DISABLED";
  this->parameters["REUSE_SOURCE"] = "false";
  this->parameters["WRITE_SOURCE"] = "false";
  // sets the kernel to verbose
  this->parameters["VERBOSE"] = "false";
  // sets the manager to verbose
  this->parameters["OCL_MANAGER_VERBOSE"] = "false";
  this->parameters["LOAD_BALANCING_VERBOSE"] = "false";
  this->parameters["SHOW_BUILD_LOG"] = "false";
}

OCLConfigurationParameters::~OCLConfigurationParameters() {}

std::shared_ptr<ConfigurationParameters> OCLConfigurationParameters::clone() {
  return std::make_shared<OCLConfigurationParameters>(*this);
  // return std::shared_ptr<ConfigurationParameters>(
  // new OCLConfigurationParameters(*this));
}

}  // namespace base
}  // namespace sgpp
