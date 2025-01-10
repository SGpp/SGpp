// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/globaldef.hpp>

#include <map>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

class OperationConfiguration : public json::JSON {
 public:
  OperationConfiguration();

  using DictNode::operator=;

  explicit OperationConfiguration(const std::string& fileName);

  virtual OperationConfiguration* clone();
};

}  // namespace base
}  // namespace sgpp
