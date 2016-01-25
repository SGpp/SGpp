/*
 * JSONIDNode.hpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/base/tools/json/Node.hpp>

namespace json {

  //enum class InternalIDType {ID, DOUBLE, UINT, INT, BOOL};

  class IDNode: public Node {
    private:
      std::string value;

      //  bool isNumber;
      //  InternalIDType internalType;
      bool isDouble;
      double doubleValue; //only used for number types
      bool isUnsigned;
      uint64_t unsignedValue;
      bool isSigned;
      int64_t signedValue;
      bool isBool;
      bool boolValue;

      void setupInternalType();

    public:
      IDNode();

      virtual void parse(std::vector<Token>& stream) override;

      virtual void serialize(std::ofstream& outFile, size_t indentWidth) override;

      virtual std::string& get() override;

      virtual void set(const std::string& value) override;

      virtual double getDouble() override;

      virtual void setDouble(double numericValue) override;

      virtual uint64_t getUInt() override;

      virtual void setUInt(uint64_t uintValue) override;

      virtual int64_t getInt() override;

      virtual void setInt(int64_t intValue) override;

      virtual bool getBool() override;

      virtual void setBool(bool boolValue) override;

      virtual size_t size() override;

      virtual Node* clone() override;
  };

}
