// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/json/DictNode.hpp>
#include <sgpp/base/tools/json/Token.hpp>

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace json {

class JSON : public DictNode {
 private:
  std::string fileName;

  std::vector<Token> tokenize(const std::string& input);

 public:
  explicit JSON(const std::string& fileName);

  JSON();

  JSON(const JSON& original);

  using DictNode::operator=;

  JSON& operator=(const JSON& original);

  virtual JSON* clone();

  void clear();

  void serialize(const std::string& outFileName);

  void deserialize(std::string content);

  void deserializeFromString(const std::string& content);

  using DictNode::serialize;
};

}  // namespace json
