// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/OperationConfiguration.hpp>

#include <string>

namespace sgpp {
namespace base {

OperationConfiguration::OperationConfiguration(): json::JSON() {
}

OperationConfiguration::OperationConfiguration(const std::string& fileName):
  json::JSON(fileName) {
}

OperationConfiguration* OperationConfiguration::clone() {
  OperationConfiguration* clone = new OperationConfiguration(*this);
  return clone;
}

}  // namespace base
}  // namespace sgpp

