// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/json/DictNode.hpp>

#include <sgpp/base/tools/json/IDNode.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>
#include <sgpp/base/tools/json/ListNode.hpp>
#include <sgpp/base/tools/json/TextNode.hpp>

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
namespace json {

DictNode::DictNode() {}

DictNode::DictNode(const DictNode& original) {
  this->keyOrder = original.keyOrder;
  this->orderedKeyIndex = original.orderedKeyIndex;
  this->parent = nullptr;

  for (auto& tuple : original.attributes) {
    std::unique_ptr<Node> clonedValue(tuple.second->clone());
    clonedValue->parent = this;
    this->attributes[tuple.first] = std::move(clonedValue);
  }
}

DictNode& DictNode::operator=(const DictNode& right) {
  this->keyOrder = right.keyOrder;
  this->orderedKeyIndex = right.orderedKeyIndex;
  this->parent = nullptr;

  this->attributes.clear();

  for (auto& tuple : right.attributes) {
    std::unique_ptr<Node> clonedValue(tuple.second->clone());
    clonedValue->parent = this;
    this->attributes[tuple.first] = std::move(clonedValue);
  }

  return *this;
}

Node& DictNode::operator=(const Node& right) {
  const DictNode& dictNode = dynamic_cast<const DictNode&>(right);
  this->operator=(dictNode);
  return *this;
}

void DictNode::parse(std::vector<Token>& stream) {
  // special case for initial and final brace
  if (stream[0].type != TokenType::LBRACE) {
    throw json_exception(stream[0], "expected \"{\"");
  }

  stream.erase(stream.begin());

  // special case for empty dict
  if (stream[0].type != TokenType::RBRACE) {
    this->parseAttributes(stream);
  }

  if (stream[0].type != TokenType::RBRACE) {
    throw json_exception(stream[0], "expected \"}\"");
  }

  stream.erase(stream.begin());
}

void DictNode::parseAttributes(std::vector<Token>& stream) {
  //  enum class Rules {
  //    NONE, TEXT_ASSIGN, ID_ASSIGN, LIST_ASSIGN, DICT_ASSIGN
  //  };
  //
  //  Rules rule = Rules::NONE;
  //  size_t i = 0;

  enum class State { NEXT = 0, COLON = 1, VALUE = 2, COMMAFINISH = 3 };

  State state = State::NEXT;

  std::string attributeName;

  while (stream.size() > 0) {
    if (state == State::NEXT) {
      if (stream[0].type == TokenType::STRING) {
        attributeName = stream[0].value;
        state = State::COLON;
        stream.erase(stream.begin());
      } else {
        throw json_exception(stream[0], "expected attribute key");
      }
    } else if (state == State::COLON) {
      if (stream[0].type == TokenType::COLON) {
        state = State::VALUE;
        stream.erase(stream.begin());
      } else {
        throw json_exception(stream[0], "expected \":\"");
      }
    } else if (state == State::VALUE) {
      if (stream[0].type == TokenType::STRING) {
        auto textNode = std::unique_ptr<TextNode>(new TextNode());
        textNode->parse(stream);
        textNode->orderedKeyIndex = this->keyOrder.size();
        this->keyOrder.push_back(attributeName);
        textNode->parent = this;
        this->attributes[attributeName] = std::move(textNode);
        state = State::COMMAFINISH;
      } else if (stream[0].type == TokenType::ID) {
        auto idNode = std::unique_ptr<IDNode>(new IDNode());
        idNode->parse(stream);
        idNode->orderedKeyIndex = this->keyOrder.size();
        this->keyOrder.push_back(attributeName);
        idNode->parent = this;
        this->attributes[attributeName] = std::move(idNode);
        state = State::COMMAFINISH;
      } else if (stream[0].type == TokenType::LBRACKET) {
        auto listNode = std::unique_ptr<ListNode>(new ListNode());
        listNode->parse(stream);
        listNode->orderedKeyIndex = this->keyOrder.size();
        this->keyOrder.push_back(attributeName);
        listNode->parent = this;
        this->attributes[attributeName] = std::move(listNode);
        state = State::COMMAFINISH;
      } else if (stream[0].type == TokenType::LBRACE) {
        auto attributeNode = std::unique_ptr<DictNode>(new DictNode());
        attributeNode->parse(stream);
        attributeNode->orderedKeyIndex = this->keyOrder.size();
        this->keyOrder.push_back(attributeName);
        attributeNode->parent = this;
        this->attributes[attributeName] = std::move(attributeNode);
        state = State::COMMAFINISH;
      } else {
        throw json_exception(stream[0],
                             "expected attribute value type "
                             "(string, id, list or dict)");
      }
    } else if (state == State::COMMAFINISH) {
      if (stream[0].type == TokenType::COMMA) {
        stream.erase(stream.begin());
        state = State::NEXT;
      } else if (stream[0].type == TokenType::RBRACE) {
        return;
      } else {
        throw json_exception(stream[0], "expected \",\" or \"}\"");
      }
    }
  }

  throw json_exception("unexpected end-of-file");
}

Node& DictNode::operator[](const std::string& key) {
  if (this->attributes.count(key) == 0) {
    throw json_exception("operator[](): key not found: " + key);
  }

  return *(this->attributes[key]);
}

void DictNode::serialize(std::ostream& outFile, size_t indentWidth) {
  std::string indentation(indentWidth, ' ');
  std::string attrIndentation(indentWidth + Node::SERIALIZE_INDENT, ' ');

  outFile << "{" << std::endl;
  bool first = true;

  for (const std::string& key : this->keyOrder) {
    if (first) {
      first = false;
    } else {
      outFile << "," << std::endl;
    }

    outFile << attrIndentation << "\"" << key << "\": ";
    this->attributes[key]->serialize(outFile, indentWidth + Node::SERIALIZE_INDENT);
  }

  outFile << std::endl << indentation << "}";
}

size_t DictNode::size() { return this->keyOrder.size(); }

Node* DictNode::clone() {
  DictNode* newNode = new DictNode(*this);
  return newNode;
}

void DictNode::addAttribute(const std::string& name, std::unique_ptr<Node> node) {
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

std::unique_ptr<Node> DictNode::removeAttribute(const std::string name) {
  if (this->attributes.count(name) == 0) {
    throw json_exception("removeAttribute(): attribute not found");
  }

  size_t orderedIndex = this->attributes[name]->orderedKeyIndex;
  auto attribute = std::move(this->attributes[name]);
  attribute->orderedKeyIndex = 0;
  attribute->parent = nullptr;
  size_t erased = this->attributes.erase(name);
  if (erased != 1) {
    throw json_exception("removeAttribute(): attribute was not erased");
  }
  this->keyOrder.erase(this->keyOrder.begin() + orderedIndex);
  // fix the ordered indices of the remaining attributes
  for (auto it = this->attributes.begin(); it != this->attributes.end(); it++) {
    Node& node = *it->second;
    if (node.orderedKeyIndex > orderedIndex) {
      node.orderedKeyIndex -= 1;
    }
  }
  return attribute;
}

// returns the node to which the attribute was added
Node& DictNode::addTextAttr(const std::string& name, const std::string& value) {
  auto textNode = std::unique_ptr<TextNode>(new TextNode());
  textNode->set(value);
  this->addAttribute(name, std::move(textNode));
  return *this;
}

// returns the node to which the attribute was added
Node& DictNode::addIDAttr(const std::string& name, const std::string& value) {
  auto idNode = std::unique_ptr<IDNode>(new IDNode());
  idNode->set(value);
  this->addAttribute(name, std::move(idNode));
  return *this;
}

// returns the node to which the attribute was added
// cast internally to string, prevents the boolean overload from being used, if the value is a
// string literal
Node& DictNode::addIDAttr(const std::string& name, const char* value) {
  return this->addIDAttr(name, std::string(value));
}

// returns the node to which the attribute was added
Node& DictNode::addIDAttr(const std::string& name, const double& value) {
  auto idNode = std::unique_ptr<IDNode>(new IDNode());
  idNode->setDouble(value);
  this->addAttribute(name, std::move(idNode));
  return *this;
}

// returns the node to which the attribute was added
Node& DictNode::addIDAttr(const std::string& name, const uint64_t& value) {
  auto idNode = std::unique_ptr<IDNode>(new IDNode());
  idNode->setUInt(value);
  this->addAttribute(name, std::move(idNode));
  return *this;
}

// returns the node to which the attribute was added
Node& DictNode::addIDAttr(const std::string& name, const int64_t& value) {
  auto idNode = std::unique_ptr<IDNode>(new IDNode());
  idNode->setInt(value);
  this->addAttribute(name, std::move(idNode));
  return *this;
}

// returns the node to which the attribute was added
Node& DictNode::addIDAttr(const std::string& name, const bool& value) {
  auto idNode = std::unique_ptr<IDNode>(new IDNode());
  idNode->setBool(value);
  this->addAttribute(name, std::move(idNode));
  return *this;
}

// returns created dict node
Node& DictNode::addDictAttr(const std::string& name) {
  auto dictNode = std::unique_ptr<DictNode>(new DictNode());
  auto& reference = *dictNode;  // because dictNode will be invalidated
  this->addAttribute(name, std::move(dictNode));
  return reference;
}

// returns created list node
Node& DictNode::addListAttr(const std::string& name) {
  auto listNode = std::unique_ptr<ListNode>(new ListNode());
  auto& reference = *listNode;  // because listNode will be invalidated
  this->addAttribute(name, std::move(listNode));
  return reference;
}

// returns the node to which the attribute was added
// replaces a node, adds a new node, if the node does not exist,
// the old node is deleted
Node& DictNode::replaceTextAttr(const std::string& name, const std::string& value) {
  if (this->attributes.count(name) > 0) {
    this->removeAttribute(name);
  }

  this->addTextAttr(name, value);
  return *this;
}

// returns the node to which the attribute was added
// replaces a node, adds a new node, if the node does not exist,
// the old node is deleted
Node& DictNode::replaceIDAttr(const std::string& name, const std::string& value) {
  if (this->attributes.count(name) > 0) {
    this->removeAttribute(name);
  }

  this->addIDAttr(name, value);
  return *this;
}

// returns the node to which the attribute was added
// replaces a node, adds a new node, if the node does not exist,
// the old node is deleted
// cast internally to string, prevents the boolean overload from being used, if the value is a
// string literal
Node& DictNode::replaceIDAttr(const std::string& name, const char* value) {
  return this->replaceIDAttr(name, std::string(value));
}

// returns the node to which the attribute was added
// replaces a node, adds a new node, if the node does not exist,
// the old node is deleted
Node& DictNode::replaceIDAttr(const std::string& name, const double& value) {
  if (this->attributes.count(name) > 0) {
    this->removeAttribute(name);
  }

  this->addIDAttr(name, value);
  return *this;
}

// returns the node to which the attribute was added
// replaces a node, adds a new node, if the node does not exist,
// the old node is deleted
Node& DictNode::replaceIDAttr(const std::string& name, const uint64_t& value) {
  if (this->attributes.count(name) > 0) {
    this->removeAttribute(name);
  }

  this->addIDAttr(name, value);
  return *this;
}

// returns the node to which the attribute was added
// replaces a node, adds a new node, if the node does not exist,
// the old node is deleted
Node& DictNode::replaceIDAttr(const std::string& name, const int64_t& value) {
  if (this->attributes.count(name) > 0) {
    this->removeAttribute(name);
  }

  this->addIDAttr(name, value);
  return *this;
}

// returns the node to which the attribute was added
// replaces a node, adds a new node, if the node does not exist,
// the old node is deleted
Node& DictNode::replaceIDAttr(const std::string& name, const bool& value) {
  if (this->attributes.count(name) > 0) {
    this->removeAttribute(name);
  }

  this->addIDAttr(name, value);
  return *this;
}

// returns created dict node
// replaces a node, adds a new node, if the node does not exist,
// the old node is deleted
Node& DictNode::replaceDictAttr(const std::string& name) {
  if (this->attributes.count(name) > 0) {
    this->removeAttribute(name);
  }

  auto& newNode = this->addDictAttr(name);
  return newNode;
}

// returns created list node
// replaces a node, adds a new node, if the node does not exist,
// the old node is deleted
Node& DictNode::replaceListAttr(const std::string& name) {
  if (this->attributes.count(name) > 0) {
    this->removeAttribute(name);
  }

  auto& newNode = this->addListAttr(name);
  return newNode;
}

bool DictNode::contains(const std::string& key) {
  std::cout<<"DictNode::contains"<<std::endl;
  std::cout<<key<<std::endl;
  if (this->attributes.count(key) > 0) {
    return true;
  } else {
    return false;
  }
}

std::unique_ptr<Node> DictNode::erase(Node& node) {
  for (auto it = this->attributes.begin(); it != this->attributes.end(); it++) {
    if (it->second.get() == &node) {
      auto temporary = this->removeAttribute(it->first);
      return temporary;
    }
  }

  throw json_exception("erase(node): node not found");
}

std::vector<std::string>& DictNode::keys() { return this->keyOrder; }

}  // namespace json
