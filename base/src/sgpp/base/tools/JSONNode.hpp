/*
 * JSONNode.hpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#pragma once

#include <vector>

#include <sgpp/globaldef.hpp>

#include "JSONToken.hpp"

namespace SGPP {
namespace base {

class JSONNode {
protected:
  static const int SERIALIZE_INDENT = 3;
  //to keep track of the insertion order for attribute nodes

public:
  size_t orderedKeyIndex;

  virtual ~JSONNode() = default;

  virtual void parse(std::vector<JSONToken> &stream) = 0;

  virtual void serialize(std::ofstream &outFile, size_t indentWidth) = 0;

  virtual JSONNode &operator[](std::string key);

  virtual std::string &getValue();

  virtual JSONNode &getItem(size_t index);

  virtual size_t size() = 0;
};

}
}
