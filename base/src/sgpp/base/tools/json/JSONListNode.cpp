/*
 * JSONListNode.cpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#include <fstream>

#include "JSONListNode.hpp"
#include "JSONIDNode.hpp"
#include "json_exception.hpp"
#include "JSONDictNode.hpp"
#include "JSONTextNode.hpp"

namespace json {

JSONListNode::JSONListNode() :
    list() {
}

void JSONListNode::parse(std::vector<JSONToken> &stream) {

  list.clear();

  enum class State {
    ITEMVALUE, NEXT
  };

  if (stream[0].type != JSONTokenType::LBRACKET) {
    throw json_exception("expected \"[\"");
  }
  stream.erase(stream.begin());

  //special case for empty list
  if (stream[0].type == JSONTokenType::RBRACKET) {
    stream.erase(stream.begin());
    return;
  }

  State state = State::ITEMVALUE;

  while (stream.size() > 0) {
    if (state == State::ITEMVALUE) {
      if (stream[0].type == JSONTokenType::STRING) {
        auto textNode = std::unique_ptr<JSONTextNode>(new JSONTextNode());
        textNode->parse(stream);
        textNode->parent = this;
        this->list.push_back(std::move(textNode));
        state = State::NEXT;
      } else if (stream[0].type == JSONTokenType::ID) {
        auto idNode = std::unique_ptr<JSONIDNode>(new JSONIDNode());
        idNode->parse(stream);
        idNode->parent = this;
        this->list.push_back(std::move(idNode));
        state = State::NEXT;
      } else if (stream[0].type == JSONTokenType::LBRACKET) {
        auto listNode = std::unique_ptr<JSONListNode>(new JSONListNode());
        listNode->parse(stream);
        listNode->parent = this;
        this->list.push_back(std::move(listNode));
        state = State::NEXT;
      } else if (stream[0].type == JSONTokenType::LBRACE) {
        auto dictNode = std::unique_ptr<JSONDictNode>(new JSONDictNode());
        dictNode->parse(stream);
        dictNode->parent = this;
        this->list.push_back(std::move(dictNode));
        state = State::NEXT;
      } else {
        throw json_exception(stream[0], "expected list value type (string, id, list or dict)");
      }
    } else if (state == State::NEXT) {
      if (stream[0].type == JSONTokenType::COMMA) {
        stream.erase(stream.begin());
        state = State::ITEMVALUE;
      } else if (stream[0].type == JSONTokenType::RBRACKET) {
        stream.erase(stream.begin());
        return;
      } else {
        throw json_exception(stream[0], "expected \",\" or \"]\"");
      }
    }
  }
  throw json_exception("unexpected end-of-file");
}

JSONNode &JSONListNode::operator[](size_t index) {
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

void JSONListNode::addValue(std::unique_ptr<JSONNode> node) {
  if (node->parent != nullptr) {
    throw json_exception("addItem(): value was already added");
  }
  node->parent = this;
  this->list.push_back(std::move(node));
}

JSONNode *JSONListNode::clone() {
  JSONListNode *newNode = new JSONListNode();
  for (auto &node : this->list) {
    newNode->list.push_back(std::unique_ptr<JSONNode>(node->clone()));
  }
  return newNode;
}

std::unique_ptr<JSONNode> JSONListNode::removeValue(size_t index) {
  if (index >= list.size()) {
    throw json_exception("removeItem(): index is out-of-bounds");
  }
  auto node = std::move(this->list[index]);
  node->parent = nullptr;
  this->list.erase(this->list.begin() + index);
  return node;
}

// returns the list node to which the value was added
JSONNode &JSONListNode::addTextValue(const std::string &value) {
  auto textNode = std::unique_ptr<JSONTextNode>(new JSONTextNode());
  textNode->set(value);
  this->addValue(std::move(textNode));
  return *this;
}

// returns the list node to which the value was added
JSONNode &JSONListNode::addIdValue(const std::string &value) {
  auto idNode = std::unique_ptr<JSONIDNode>(new JSONIDNode());
  idNode->set(value);
  this->addValue(std::move(idNode));
  return *this;
}

// returns the list node to which the value was added
JSONNode &JSONListNode::addIdValue(const double &value) {
  auto idNode = std::unique_ptr<JSONIDNode>(new JSONIDNode());
  idNode->setNumeric(value);
  this->addValue(std::move(idNode));
  return *this;
}

// returns created dict node
JSONNode &JSONListNode::addDictValue() {
  auto dictNode = std::unique_ptr<JSONDictNode>(new JSONDictNode());
  auto &reference = *dictNode; // because dictNode will be invalidated
  this->addValue(std::move(dictNode));
  return reference;
}

// returns created dict node
JSONNode &JSONListNode::addListValue() {
  auto listNode = std::unique_ptr<JSONListNode>(new JSONListNode());
  auto &reference = *listNode; // because listNode will be invalidated
  this->addValue(std::move(listNode));
  return reference;
}

std::unique_ptr<JSONNode> JSONListNode::erase(JSONNode &node) {
  for (auto it = this->list.begin(); it != this->list.end(); it++) {
    if ((*it).get() == &node) {
      auto node = std::move(*it);
      node->parent = nullptr;
      this->list.erase(it);
      return node;
    }
  }

  throw json_exception("erase(node): node not found");
}

}
