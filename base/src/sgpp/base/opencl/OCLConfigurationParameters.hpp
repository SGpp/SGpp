// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <vector>
#include <map>
#include <string>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/tools/ConfigurationParameters.hpp>

namespace sgpp {
namespace base {

class OCLConfigurationParameters : public ConfigurationParameters {
 public:
  OCLConfigurationParameters(
      std::string fileName,
      std::map<std::string, std::string> defaultParameters);

  OCLConfigurationParameters();

  virtual ~OCLConfigurationParameters();

  std::shared_ptr<ConfigurationParameters> clone() override;
};

}  // namespace base
}  // namespace sgpp
