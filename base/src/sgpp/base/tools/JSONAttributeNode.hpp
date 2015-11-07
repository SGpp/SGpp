/*
 * JSONAttributeNode.hpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#pragma once

#include <map>
#include <memory>

#include <sgpp/globaldef.hpp>

#include "JSONNode.hpp"

namespace SGPP {
namespace base {

class JSONAttributeNode: public JSONNode {
private:
  std::map<std::string, std::unique_ptr<JSONNode>> attributes;

  std::vector<std::string> keyOrder;
public:
  virtual void parse(std::vector<JSONToken> &stream) override;

  void parseAttributes(std::vector<JSONToken> &stream);

  virtual void serialize(std::ofstream &outFile, size_t indentWidth) override;

  virtual JSONNode &operator[](std::string key) override;

  virtual size_t size() override;
};

}
}
