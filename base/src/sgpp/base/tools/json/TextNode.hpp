// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/json/Node.hpp>

#include <string>
#include <vector>

namespace json {

class TextNode: public Node {
 private:
  std::string value;

 public:
  TextNode();

  void parse(std::vector<Token>& stream) override;

  void serialize(std::ofstream& outFile, size_t indentWidth) override;

  std::string& get() override;

  void set(const std::string& value) override;

  size_t size() override;

  Node* clone() override;
};

}  // namespace json
