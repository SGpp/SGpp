/*
 * JSONToken.hpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#pragma once

namespace SGPP {
namespace base {

enum class JSONTokenType {
  LBRACE, RBRACE, LBRACKET, RBRACKET, COLON, ID, STRING, COMMA, NONE
};

class JSONToken {
public:
  JSONTokenType type;
  std::string value;
};

}
}
