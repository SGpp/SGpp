// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/json/JSON.hpp>

#include <vector>
#include <map>
#include <string>

namespace SGPP {
namespace base {

class OperationConfiguration: public json::JSON {
 public:
  OperationConfiguration();

  explicit OperationConfiguration(const std::string& fileName);

  virtual OperationConfiguration* clone();
};

}  // namespace base
}  // namespace SGPP

