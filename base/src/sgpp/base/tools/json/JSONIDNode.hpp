/*
 * JSONIDNode.hpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#pragma once

#include "JSONNode.hpp"

namespace json {

class JSONIDNode: public JSONNode {
private:
  std::string value;

  bool isNumber;
  double numericValue; //only used for number types

  void tryInterpretAsNumber();

public:
  JSONIDNode();

  virtual void parse(std::vector<JSONToken> &stream) override;

  virtual void serialize(std::ofstream &outFile, size_t indentWidth) override;

  virtual std::string &get() override;

  virtual void set(const std::string &value) override;

  virtual double getNumeric() override;

  virtual void setNumeric(double numericValue) override;

  virtual size_t size() override;

  virtual JSONNode *clone() override;
};

}
