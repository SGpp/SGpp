// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/tools/ConfigurationParameters.hpp>

#include <vector>
#include <map>
#include <string>

namespace SGPP {
namespace base {

class OpenCLConfigurationParameters: public ConfigurationParameters {
 public:
  OpenCLConfigurationParameters(
    std::string fileName,
    std::map<std::string, std::string> defaultParameters);
};

}  // namespace base
}  // namespace SGPP

