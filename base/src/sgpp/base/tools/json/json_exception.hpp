// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/json/Token.hpp>

#include <exception>
#include <string>

namespace json {

class json_exception: public std::exception {
 private:
  std::string message;

 public:
  json_exception(Token& token, const std::string& message);

  explicit json_exception(const std::string& message);

  const char* what() const noexcept override;
};

}  // namespace json
