// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/json/Node.hpp>

#include <string>
#include <vector>

namespace json {

// enum class InternalIDType {ID, DOUBLE, UINT, INT, BOOL};

class IDNode: public Node {
 private:
  std::string value;

  //  bool isNumber;
  //  InternalIDType internalType;
  bool isDouble;
  double doubleValue;  // only used for number types
  bool isUnsigned;
  uint64_t unsignedValue;
  bool isSigned;
  int64_t signedValue;
  bool isBool;
  bool boolValue;

  void setupInternalType();

 public:
  IDNode();

  IDNode& operator=(const IDNode& right) = default;

  Node& operator=(const Node& right) override;

  void parse(std::vector<Token>& stream) override;

  void serialize(std::ostream& outFile, size_t indentWidth) override;

  std::string& get() override;

  void set(const std::string& value) override;

  double getDouble() override;

  void setDouble(double numericValue) override;

  uint64_t getUInt() override;

  void setUInt(uint64_t uintValue) override;

  int64_t getInt() override;

  void setInt(int64_t intValue) override;

  bool getBool() override;

  void setBool(bool boolValue) override;

  size_t size() override;

  Node* clone() override;
};

}  // namespace json
