/*
 * JSONListNode.hpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#pragma once

#include <vector>

#include <sgpp/globaldef.hpp>

#include "JSONNode.hpp"

namespace SGPP {
namespace base {

class JSONListNode: public JSONNode {
private:
  std::vector<std::unique_ptr<JSONNode>> list;
public:
  JSONListNode();

  void parse(std::vector<JSONToken> &stream) override;

  virtual void serialize(std::ofstream &outFile, size_t indentWidth) override;

  virtual JSONNode &getItem(size_t index) override;

  virtual size_t size() override;
};

}
}


