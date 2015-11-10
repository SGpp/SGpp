/*
 * JSONNode.hpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#pragma once

#include <memory>
#include <vector>

#include "JSONToken.hpp"

namespace json {

class JSONNode {
protected:

  static const int SERIALIZE_INDENT = 3;

public:
  size_t orderedKeyIndex;

  JSONNode *parent;

  virtual ~JSONNode() = default;

  virtual void parse(std::vector<JSONToken> &stream) = 0;

  virtual void serialize(std::ofstream &outFile, size_t indentWidth) = 0;

  virtual JSONNode &operator[](std::string key);

  virtual JSONNode &operator[](size_t index);

  virtual std::string &get();

  virtual void set(const std::string &value);

  virtual double getNumeric();

  virtual void setNumeric(double numericValue);

//  virtual JSONNode &getItem(size_t index);

  virtual size_t size() = 0;

  virtual void addValue(std::unique_ptr<JSONNode> node);

  virtual void addAttribute(const std::string &name, std::unique_ptr<JSONNode> node);

  virtual std::unique_ptr<JSONNode> removeValue(size_t index);

  virtual std::unique_ptr<JSONNode> removeAttribute(const std::string &name);

  virtual JSONNode *clone() = 0;

  // returns the node to which the attribute was added
  virtual JSONNode &addTextAttr(const std::string &name, const std::string &value);

  // returns the node to which the attribute was added
  virtual JSONNode &addIDAttr(const std::string &name, const std::string &value);

  // returns the node to which the attribute was added
  virtual JSONNode &addIDAttr(const std::string &name, const double &numericValue);

  // returns created dict node
  virtual JSONNode &addDictAttr(const std::string &name);

  // returns created list node
  virtual JSONNode &addListAttr(const std::string &name);

  // returns created dict node
  virtual JSONNode &addDictValue();

  // returns created dict node
  virtual JSONNode &addListValue();

  // returns the list node to which the value was added
  virtual JSONNode &addTextValue(const std::string &value);

  // returns the list node to which the value was added
  virtual JSONNode &addIdValue(const std::string &value);

  // returns the list node to which the value was added
  virtual JSONNode &addIdValue(const double &value);

  virtual std::unique_ptr<JSONNode> erase(JSONNode &node);

  virtual std::unique_ptr<JSONNode> erase();

};

}
