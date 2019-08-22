// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/json/json_exception.hpp>

#include <string>
#include <sstream>

namespace json {

json_exception::json_exception(Token& token, const std::string& message) {
  std::stringstream messageStream;
  messageStream << "error: (line: " << token.lineNumber << ", char: " <<
                token.charNumber << ") at \"" << token.value
                << "\": ";
  messageStream << message;
  this->message = messageStream.str();
}

json_exception::json_exception(const std::string& message) {
  std::stringstream messageStream;
  messageStream << "error: ";
  messageStream << message;
  this->message = messageStream.str();
}

const char* json_exception::what() const noexcept {
  return this->message.c_str();
}

}  // namespace json
