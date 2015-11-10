/*
 * JSONDictNode.cpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#include "JSONDictNode.hpp"

#include <fstream>

#include "JSONIDNode.hpp"
#include "JSONListNode.hpp"
#include "json_exception.hpp"
#include "JSONTextNode.hpp"

namespace json {

void JSONDictNode::parse(std::vector<JSONToken> &stream) {

  //special case for initial and final brace
  if (stream[0].type != JSONTokenType::LBRACE) {
    throw json_exception(stream[0], "expected \"{\"");
  }
  stream.erase(stream.begin());

  this->parseAttributes(stream);

  if (stream[0].type != JSONTokenType::RBRACE) {
    throw json_exception(stream[0], "expected \"}\"");
  }
  stream.erase(stream.begin());
}

void JSONDictNode::parseAttributes(std::vector<JSONToken> &stream) {
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
        throw json_exception(stream[0], "expected attribute key");
      }
    } else if (state == State::COLON) {
      if (stream[0].type == JSONTokenType::COLON) {
        state = State::VALUE;
        stream.erase(stream.begin());
      } else {
        throw json_exception(stream[0], "expected \":\"");
      }
    } else if (state == State::VALUE) {
      if (stream[0].type == JSONTokenType::STRING) {
        auto textNode = std::unique_ptr<JSONTextNode>(new JSONTextNode());
        textNode->parse(stream);
        textNode->orderedKeyIndex = this->keyOrder.size();
        this->keyOrder.push_back(attributeName);
        textNode->parent = this;
        this->attributes[attributeName] = std::move(textNode);
        state = State::COMMAFINISH;
      } else if (stream[0].type == JSONTokenType::ID) {
        auto idNode = std::unique_ptr<JSONIDNode>(new JSONIDNode());
        idNode->parse(stream);
        idNode->orderedKeyIndex = this->keyOrder.size();
        this->keyOrder.push_back(attributeName);
        idNode->parent = this;
        this->attributes[attributeName] = std::move(idNode);
        state = State::COMMAFINISH;
      } else if (stream[0].type == JSONTokenType::LBRACKET) {
        auto listNode = std::unique_ptr<JSONListNode>(new JSONListNode());
        listNode->parse(stream);
        listNode->orderedKeyIndex = this->keyOrder.size();
        this->keyOrder.push_back(attributeName);
        listNode->parent = this;
        this->attributes[attributeName] = std::move(listNode);
        state = State::COMMAFINISH;
      } else if (stream[0].type == JSONTokenType::LBRACE) {
        auto attributeNode = std::unique_ptr<JSONDictNode>(new JSONDictNode());
        attributeNode->parse(stream);
        attributeNode->orderedKeyIndex = this->keyOrder.size();
        this->keyOrder.push_back(attributeName);
        attributeNode->parent = this;
        this->attributes[attributeName] = std::move(attributeNode);
        state = State::COMMAFINISH;
      } else {
        throw json_exception(stream[0], "expected attribute value type (string, id, list or dict)");
      }
    } else if (state == State::COMMAFINISH) {
      if (stream[0].type == JSONTokenType::COMMA) {
        stream.erase(stream.begin());
        state = State::NEXT;
      } else if (stream[0].type == JSONTokenType::RBRACE) {
        return;
      } else {
        throw json_exception(stream[0], "expected \",\" or \"}\"");
      }
    }
  }

  throw json_exception("unexpected end-of-file");

}

JSONNode &JSONDictNode::operator[](std::string key) {
  return *(this->attributes[key]);
}

void JSONDictNode::serialize(std::ofstream &outFile, size_t indentWidth) {
  std::string indentation(indentWidth, ' ');
  std::string attrIndentation(indentWidth + JSONNode::SERIALIZE_INDENT, ' ');

  outFile << "{" << std::endl;
  bool first = true;
  for (std::string &key : this->keyOrder) {
    if (first) {
      first = false;
    } else {
      outFile << "," << std::endl;
    }
    outFile << attrIndentation << "\"" << key << "\": ";
    this->attributes[key]->serialize(outFile, indentWidth + JSONNode::SERIALIZE_INDENT);

  }

  outFile << std::endl << indentation << "}";
}

size_t JSONDictNode::size() {
  return this->keyOrder.size();
}

JSONNode *JSONDictNode::clone() {

  JSONDictNode *newNode = new JSONDictNode();

  for (auto &tuple: this->attributes) {
    newNode->attributes[tuple.first] = std::unique_ptr<JSONNode>(tuple.second->clone());
  }
  newNode->keyOrder = this->keyOrder;

  return newNode;
}

void JSONDictNode::addAttribute(const std::string &name, std::unique_ptr<JSONNode> node) {
  if (node->parent != nullptr) {
    throw json_exception("addAttribute(): attribute was already added");
  } else if (this->attributes.count(name) > 0) {
    throw json_exception("addAttribute(): attribute with same name already exists");
  }
  node->parent = this;
  node->orderedKeyIndex = this->keyOrder.size();
  this->attributes[name] = std::move(node);
  this->keyOrder.push_back(name);
}

std::unique_ptr<JSONNode> JSONDictNode::removeAttribute(const std::string &name) {
  if (this->attributes.count(name) == 0) {
    throw json_exception("removeAttribute(): attribute not found");
  }
  size_t orderedIndex = this->attributes[name]->orderedKeyIndex;
  auto attribute = std::move(this->attributes[name]);
  attribute->orderedKeyIndex = 0;
  attribute->parent = nullptr;
  this->attributes.erase(name);
  this->keyOrder.erase(this->keyOrder.begin() + orderedIndex);
  return attribute;
}

// returns the node to which the attribute was added
JSONNode &JSONDictNode::addTextAttr(const std::string &name, const std::string &value) {
  auto textNode = std::unique_ptr<JSONTextNode>(new JSONTextNode());
  textNode->set(value);
  this->addAttribute(name, std::move(textNode));
  return *this;
}

// returns the node to which the attribute was added
JSONNode &JSONDictNode::addIDAttr(const std::string &name, const std::string &value) {
  auto idNode = std::unique_ptr<JSONIDNode>(new JSONIDNode());
  idNode->set(value);
  this->addAttribute(name, std::move(idNode));
  return *this;
}

// returns the node to which the attribute was added
JSONNode &JSONDictNode::addIDAttr(const std::string &name, const double &numericValue) {
  auto idNode = std::unique_ptr<JSONIDNode>(new JSONIDNode());
  idNode->setNumeric(numericValue);
  this->addAttribute(name, std::move(idNode));
  return *this;
}

// returns created dict node
JSONNode &JSONDictNode::addDictAttr(const std::string &name) {
  auto dictNode = std::unique_ptr<JSONDictNode>(new JSONDictNode());
  auto &reference = *dictNode; // because dictNode will be invalidated
  this->addAttribute(name, std::move(dictNode));
  return reference;
}

// returns created list node
JSONNode &JSONDictNode::addListAttr(const std::string &name) {
  auto listNode = std::unique_ptr<JSONListNode>(new JSONListNode());
  auto &reference = *listNode; // because listNode will be invalidated
  this->addAttribute(name, std::move(listNode));
  return reference;
}

std::unique_ptr<JSONNode> JSONDictNode::erase(JSONNode &node) {
  for (auto it = this->attributes.begin(); it != this->attributes.end(); it++) {
    if (it->second.get() == &node) {
      auto temporary = this->removeAttribute(it->first);
      return temporary;
    }
  }

  throw json_exception("erase(node): node not found");
}

}
