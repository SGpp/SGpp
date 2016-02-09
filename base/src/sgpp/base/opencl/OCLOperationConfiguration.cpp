// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

#include <string>

namespace SGPP {
namespace base {

OCLOperationConfiguration::OCLOperationConfiguration():
  OperationConfiguration() {
  // augment default values to configuration
  //        this->addIDAttr("LOCAL_SIZE", 64ul);
  //        this->addIDAttr("ENABLE_OPTIMIZATIONS", true);
  //        this->addTextAttr("OPTIMIZATION_FLAGS", "");
  //        this->addTextAttr("INTERNAL_PRECISION", "double");
  //        this->addTextAttr("PLATFORM", "first");
  //        this->addTextAttr("DEVICE_TYPE", "CL_DEVICE_TYPE_ALL");
  //        this->addIDAttr("REUSE_SOURCE", false);
  //        this->addIDAttr("WRITE_SOURCE", false);
  //        // sets the kernel to verbose
  //        this->addIDAttr("VERBOSE", true);
  //        // sets the manager to verbose
  //        this->addIDAttr("OCL_MANAGER_VERBOSE", false);
  ////        this->addIDAttr("LOAD_BALANCING_VERBOSE", false);
  //        this->addIDAttr("SHOW_BUILD_LOG", false);
}

OCLOperationConfiguration::OCLOperationConfiguration(const std::string&
    fileName): OperationConfiguration(fileName) {
  // augment default values to configuration
  //    if (this->contains("LOCAL_SIZE") == false) {
  //        this->addIDAttr("LOCAL_SIZE", 64ul);
  //    }
  //
  //    if (this->contains("ENABLE_OPTIMIZATIONS") == false) {
  //        this->addIDAttr("ENABLE_OPTIMIZATIONS", true);
  //    }
  //
  //    if (this->contains("OPTIMIZATION_FLAGS") == false) {
  //        this->addTextAttr("OPTIMIZATION_FLAGS", "");
  //    }
  //
  //    if (this->contains("INTERNAL_PRECISION") == false) {
  //        this->addTextAttr("INTERNAL_PRECISION", "double");
  //    }
  //
  //    if (this->contains("PLATFORM") == false) {
  //        this->addTextAttr("PLATFORM", "first");
  //    }
  //
  //    if (this->contains("DEVICE_TYPE") == false) {
  //        this->addTextAttr("DEVICE_TYPE", "CL_DEVICE_TYPE_ALL");
  //    }
  //
  //    if (this->contains("REUSE_SOURCE") == false) {
  //        this->addIDAttr("REUSE_SOURCE", false);
  //    }
  //
  //    if (this->contains("WRITE_SOURCE") == false) {
  //        this->addIDAttr("WRITE_SOURCE", false);
  //    }
  //
  //    if (this->contains("VERBOSE") == false) {
  //        // sets the kernel to verbose
  //        this->addIDAttr("VERBOSE", true);
  //    }
  //
  //    if (this->contains("OCL_MANAGER_VERBOSE") == false) {
  //        // sets the manager to verbose
  //        this->addIDAttr("OCL_MANAGER_VERBOSE", false);
  //    }
  //
  //    if (this->contains("LOAD_BALANCING_VERBOSE") == false) {
  //        this->addIDAttr("LOAD_BALANCING_VERBOSE", false);
  //    }
  //
  //    if (this->contains("SHOW_BUILD_LOG") == false) {
  //        this->addIDAttr("SHOW_BUILD_LOG", false);
  //    }
}

OperationConfiguration* OCLOperationConfiguration::clone() {
  return dynamic_cast<OperationConfiguration*>(new OCLOperationConfiguration(
           *this));
}


}  // namespace base
}  // namespace SGPP

