// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/OperationConfiguration.hpp>

#include <vector>
#include <map>
#include <string>

namespace SGPP {
namespace base {

class OCLOperationConfiguration: public OperationConfiguration {
 public:
  OCLOperationConfiguration();

  explicit OCLOperationConfiguration(const std::string& fileName);

  OperationConfiguration* clone() override;
};

}  // namespace base
}  // namespace SGPP

