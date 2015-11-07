/*
 * JSONNode.cpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#include "JSONNode.hpp"

namespace SGPP {
namespace base {


JSONNode &JSONNode::operator[](std::string key) {
  throw;
}

std::string &JSONNode::getValue() {
  throw;
}

JSONNode &JSONNode::getItem(size_t index) {
  throw;
}

}
}
