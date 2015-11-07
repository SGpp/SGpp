/*
 * JSONListNode.cpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#include <fstream>

#include "JSONListNode.hpp"
#include "JSONStringNode.hpp"
#include "JSONIDNode.hpp"
#include "JSONAttributeNode.hpp"

namespace SGPP {
namespace base {

JSONListNode::JSONListNode() :
    list() {
}

void JSONListNode::parse(std::vector<JSONToken> &stream) {

  list.clear();

  enum class State {
    ITEMVALUE, NEXT
  };

  if (stream[0].type != JSONTokenType::LBRACKET) {
    throw;
  }
  stream.erase(stream.begin());

  //special case for empty list
  if (stream[0].type == JSONTokenType::RBRACKET) {
    stream.erase(stream.begin());
    return;
  }

  State state = State::ITEMVALUE;

  while (true) {
    if (state == State::ITEMVALUE) {
      if (stream[0].type == JSONTokenType::STRING) {
        auto textNode = std::make_unique<JSONStringNode>();
        textNode->parse(stream);
        this->list.push_back(std::move(textNode));
        state = State::NEXT;
      } else if (stream[0].type == JSONTokenType::ID) {
        auto idNode = std::make_unique<JSONIDNode>();
        idNode->parse(stream);
        this->list.push_back(std::move(idNode));
        state = State::NEXT;
      } else if (stream[0].type == JSONTokenType::LBRACKET) {
        auto listNode = std::make_unique<JSONListNode>();
        listNode->parse(stream);
        this->list.push_back(std::move(listNode));
        state = State::NEXT;
      } else if (stream[0].type == JSONTokenType::LBRACE) {
        auto attributeNode = std::make_unique<JSONAttributeNode>();
        attributeNode->parse(stream);
        this->list.push_back(std::move(attributeNode));
        state = State::NEXT;
      } else {
        throw; //expected a value type
      }

//      if (stream[0].type == JSONTokenType::STRING) {
//        list.push_back(stream[0].value);
//        stream.erase(stream.begin());
//        state = State::NEXT;
//      } else {
//        throw;
//      }
    } else if (state == State::NEXT) {
      if (stream[0].type == JSONTokenType::COMMA) {
        stream.erase(stream.begin());
        state = State::ITEMVALUE;
      } else if (stream[0].type == JSONTokenType::RBRACKET) {
        stream.erase(stream.begin());
        return;
      } else {
        throw;
      }
    }
  }
  throw;
}

JSONNode &JSONListNode::getItem(size_t index) {
  return *this->list[index];
}

void JSONListNode::serialize(std::ofstream &outFile, size_t indentWidth) {
  outFile << "[";
  bool first = true;
  for (const std::unique_ptr<JSONNode> &node : this->list) {
    if (first) {
      first = false;
    } else {
      outFile << ", ";
    }
    node->serialize(outFile, indentWidth + JSONNode::SERIALIZE_INDENT);
  }
  outFile << "]";
}

size_t JSONListNode::size() {
  return this->list.size();
}

}
}

