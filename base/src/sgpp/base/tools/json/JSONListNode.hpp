/*
 * JSONListNode.hpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#pragma once

#include <memory>
#include <vector>

#include "JSONNode.hpp"

namespace json {

class JSONListNode: public JSONNode {
private:

  std::vector<std::unique_ptr<JSONNode>> list;

public:
  JSONListNode();

  void parse(std::vector<JSONToken> &stream) override;

  virtual void serialize(std::ofstream &outFile, size_t indentWidth) override;

  virtual JSONNode &operator[](size_t index) override;

  virtual size_t size() override;

  virtual void addValue(std::unique_ptr<JSONNode> node) override;

  virtual std::unique_ptr<JSONNode> removeValue(size_t index) override;

  virtual JSONNode *clone() override;

  // returns created dict node
  virtual JSONNode &addDictValue() override;

  // returns created dict node
  virtual JSONNode &addListValue() override;

  // returns the list node to which the value was added
  virtual JSONNode &addTextValue(const std::string &value) override;

  // returns the list node to which the value was added
  virtual JSONNode &addIdValue(const std::string &value) override;

  // returns the list node to which the value was added
  virtual JSONNode &addIdValue(const double &value) override;

  virtual std::unique_ptr<JSONNode> erase(JSONNode &node) override;

};

}
