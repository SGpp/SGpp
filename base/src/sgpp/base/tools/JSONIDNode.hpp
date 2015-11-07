/*
 * JSONIDNode.hpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/globaldef.hpp>

#include "JSONNode.hpp"

namespace SGPP {
namespace base {

class JSONIDNode: public JSONNode {
private:
  std::string value;
public:
  JSONIDNode();

  virtual void parse(std::vector<JSONToken> &stream) override;

  virtual void serialize(std::ofstream &outFile, size_t indentWidth) override;

  virtual std::string &getValue() override;

  virtual size_t size() override;
};

}
}
