/*
 * JSONAttributeNode.cpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#include <fstream>

#include "JSONAttributeNode.hpp"
#include "JSONIDNode.hpp"
#include "JSONStringNode.hpp"
#include "JSONListNode.hpp"

namespace SGPP {
namespace base {

void JSONAttributeNode::parse(std::vector<JSONToken> &stream) {

  //special case for initial and final brace
  if (stream[0].type != JSONTokenType::LBRACE) {
    throw;
  }
  stream.erase(stream.begin());

  this->parseAttributes(stream);

  if (stream[0].type != JSONTokenType::RBRACE) {
    throw;
  }
  stream.erase(stream.begin());
}

void JSONAttributeNode::parseAttributes(std::vector<JSONToken> &stream) {
//  enum class Rules {
//    NONE, TEXT_ASSIGN, ID_ASSIGN, LIST_ASSIGN, DICT_ASSIGN
//  };
//
//  Rules rule = Rules::NONE;
//  size_t i = 0;

  enum class State {
    NEXT = 0, COLON = 1, VALUE = 2, COMMAFINISH = 3
  };

  State state = State::NEXT;

  std::string attributeName;

  while (stream.size() > 0) {
    if (state == State::NEXT) {
      if (stream[0].type == JSONTokenType::STRING) {
        attributeName = stream[0].value;
        state = State::COLON;
        stream.erase(stream.begin());
      } else {
        throw; //expected string or end of attribute
      }
    } else if (state == State::COLON) {
      if (stream[0].type == JSONTokenType::COLON) {
        state = State::VALUE;
        stream.erase(stream.begin());
      } else {
        throw; //expected colon
      }
    } else if (state == State::VALUE) {
      if (stream[0].type == JSONTokenType::STRING) {
        auto textNode = std::make_unique<JSONStringNode>();
        textNode->parse(stream);
        textNode->orderedKeyIndex = this->keyOrder.size();
        this->keyOrder.push_back(attributeName);
        this->attributes[attributeName] = std::move(textNode);
        state = State::COMMAFINISH;
      } else if (stream[0].type == JSONTokenType::ID) {
        auto textNode = std::make_unique<JSONIDNode>();
        textNode->parse(stream);
        textNode->orderedKeyIndex = this->keyOrder.size();
        this->keyOrder.push_back(attributeName);
        this->attributes[attributeName] = std::move(textNode);
        state = State::COMMAFINISH;
      } else if (stream[0].type == JSONTokenType::LBRACKET) {
        auto listNode = std::make_unique<JSONListNode>();
        listNode->parse(stream);
        listNode->orderedKeyIndex = this->keyOrder.size();
        this->keyOrder.push_back(attributeName);
        this->attributes[attributeName] = std::move(listNode);
        state = State::COMMAFINISH;
      } else if (stream[0].type == JSONTokenType::LBRACE) {
        auto attributeNode = std::make_unique<JSONAttributeNode>();
        attributeNode->parse(stream);
        attributeNode->orderedKeyIndex = this->keyOrder.size();
        this->keyOrder.push_back(attributeName);
        this->attributes[attributeName] = std::move(attributeNode);
        state = State::COMMAFINISH;
      } else {
        throw; //expected colon
      }
    } else if (state == State::COMMAFINISH) {
      if (stream[0].type == JSONTokenType::COMMA) {
        stream.erase(stream.begin());
        state = State::NEXT;
      } else if (stream[0].type == JSONTokenType::RBRACE) {
        return;
      } else {
        throw; //expected comma or rbrace
      }
    }
  }

  throw;

}

JSONNode &JSONAttributeNode::operator[](std::string key) {
//  JSONAttributeNode *casted = const_cast<JSONAttributeNode *>(this);
//  casted->attributes[key];
//  JSONNode &subNode = *(const_cast<JSONAttributeNode *>(this)->attributes[key]);
  return *(this->attributes[key]);
}

void JSONAttributeNode::serialize(std::ofstream &outFile, size_t indentWidth) {
  std::string indentation(indentWidth, ' ');
  std::string attrIndentation(indentWidth + JSONNode::SERIALIZE_INDENT, ' ');

  outFile << indentation << "{" << std::endl;
  bool first = true;
  for (std::string &key : this->keyOrder) {
    if (first) {
      first = false;
    } else {
      outFile << "," << std::endl;
    }
    outFile << attrIndentation << "\"" << key << ": ";
    this->attributes[key]->serialize(outFile, indentWidth + JSONNode::SERIALIZE_INDENT);

  }

  outFile << std::endl << indentation << "}";
}

size_t JSONAttributeNode::size() {
  return this->keyOrder.size();
}

}
}
