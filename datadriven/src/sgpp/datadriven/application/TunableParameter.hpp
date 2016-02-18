// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <string>
#include <vector>

#include "sgpp/globaldef.hpp"

namespace SGPP {
namespace datadriven {

// enum class ParameterType {
//    ID, TEXT
//};

class TunableParameter {
 private:
  std::string name;
  std::vector<std::string> values;
  //    ParameterType type;
 public:
  //    TunableParameter(std::string name, std::vector<std::string> values, ParameterType type) :
  TunableParameter(std::string name, std::vector<std::string> values)
      :  //            name(name), values(values), type(type) {
        name(name),
        values(values) {}

  std::string &getName() { return this->name; }

  std::vector<std::string> &getValues() { return this->values; }
};
}  // namespace datadriven
}  // namespace SGPP
