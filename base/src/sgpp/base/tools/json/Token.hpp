/*
 * JSONToken.hpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#pragma once

#include <string>

namespace json {

  enum class TokenType {
    LBRACE, RBRACE, LBRACKET, RBRACKET, COLON, ID, STRING, COMMA, COMMENT, SINGLELINE, MULTILINECOMMENT, MULTILINECOMMENTSTAR, NONE
  };

  class Token {
    public:
      TokenType type;
      std::string value;
      size_t lineNumber;
      size_t charNumber;
  };

}
