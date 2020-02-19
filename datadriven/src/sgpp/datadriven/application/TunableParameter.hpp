// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

class TunableParameter {
 private:
  std::string name;
  std::vector<std::string> values;

 public:
  TunableParameter(std::string name, std::vector<std::string> values)
      : name(name), values(values) {}

  std::string &getName() { return this->name; }

  std::vector<std::string> &getValues() { return this->values; }
};
}  // namespace datadriven
}  // namespace sgpp
