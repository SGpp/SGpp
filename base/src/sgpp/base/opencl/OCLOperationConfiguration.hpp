// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <vector>
#include <map>
#include <string>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/tools/OperationConfiguration.hpp>

#if USE_OCL == 0
#error "OCL is not enabled"
#endif
#ifndef USE_OCL
#error "OCL is not defined"
#endif

namespace sgpp {
namespace base {

class OCLOperationConfiguration : public OperationConfiguration {
 public:
  OCLOperationConfiguration();

  explicit OCLOperationConfiguration(const std::string& fileName);

  OperationConfiguration* clone() override;

  std::vector<std::reference_wrapper<json::Node>> getAllDeviceNodes();

  static std::unique_ptr<OCLOperationConfiguration> fromString(std::string &parametersString);
};

}  // namespace base
}  // namespace sgpp
