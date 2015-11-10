/*
 * JSONDictNode.hpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#pragma once

#include <map>
#include <memory>

#include "JSONNode.hpp"

namespace json {

class JSONDictNode: public JSONNode {
private:
  std::map<std::string, std::unique_ptr<JSONNode>> attributes;

  std::vector<std::string> keyOrder;

public:

  virtual void parse(std::vector<JSONToken> &stream) override;

  void parseAttributes(std::vector<JSONToken> &stream);

  virtual void serialize(std::ofstream &outFile, size_t indentWidth) override;

  virtual JSONNode &operator[](std::string key) override;

  virtual size_t size() override;

  virtual JSONNode *clone() override;

  virtual void addAttribute(const std::string &name, std::unique_ptr<JSONNode> node) override;

  virtual std::unique_ptr<JSONNode> removeAttribute(const std::string &name) override;

  // returns the node to which the attribute was added
  virtual JSONNode &addTextAttr(const std::string &name, const std::string &value) override;

  // returns the node to which the attribute was added
  virtual JSONNode &addIDAttr(const std::string &name, const std::string &value) override;

  // returns the node to which the attribute was added
  virtual JSONNode &addIDAttr(const std::string &name, const double &numericValue) override;

  // returns created dict node
  virtual JSONNode &addDictAttr(const std::string &name) override;

  // returns created list node
  virtual JSONNode &addListAttr(const std::string &name) override;

  virtual std::unique_ptr<JSONNode> erase(JSONNode &node) override;
};

}
