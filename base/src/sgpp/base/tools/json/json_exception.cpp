/*
 * JSONException.hpp
 *
 *  Created on: Nov 9, 2015
 *      Author: pfandedd
 */

#include <string>
#include <sstream>

#include "json_exception.hpp"

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

const char* json_exception::what() const throw () {
  return this->message.c_str();
}

}
