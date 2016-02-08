// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/json/ListNode.hpp>

#include <sgpp/base/tools/json/DictNode.hpp>
#include <sgpp/base/tools/json/IDNode.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>
#include <sgpp/base/tools/json/TextNode.hpp>

#include <fstream>
#include <string>
#include <vector>

namespace json {

ListNode::ListNode() :
  list() {
}

ListNode::ListNode(const ListNode& original) {
  for (auto& element : original.list) {
    std::unique_ptr<Node> cloned(element->clone());
    cloned->parent = this;
    this->list.push_back(std::move(cloned));
  }

  this->orderedKeyIndex = original.orderedKeyIndex;
  this->parent = nullptr;
}

ListNode& ListNode::operator=(const ListNode& right) {
  for (auto& element : right.list) {
    std::unique_ptr<Node> cloned(element->clone());
    cloned->parent = this;
    this->list.push_back(std::move(cloned));
  }

  this->orderedKeyIndex = right.orderedKeyIndex;
  this->parent = nullptr;
  return *this;
}

void ListNode::parse(std::vector<Token>& stream) {
  list.clear();

  enum class State {
    ITEMVALUE, NEXT
  };

  if (stream[0].type != TokenType::LBRACKET) {
    throw json_exception("expected \"[\"");
  }

  stream.erase(stream.begin());

  // special case for empty list
  if (stream[0].type == TokenType::RBRACKET) {
    stream.erase(stream.begin());
    return;
  }

  State state = State::ITEMVALUE;

  while (stream.size() > 0) {
    if (state == State::ITEMVALUE) {
      if (stream[0].type == TokenType::STRING) {
        auto textNode = std::unique_ptr<TextNode>(new TextNode());
        textNode->parse(stream);
        textNode->parent = this;
        this->list.push_back(std::move(textNode));
        state = State::NEXT;
      } else if (stream[0].type == TokenType::ID) {
        auto idNode = std::unique_ptr<IDNode>(new IDNode());
        idNode->parse(stream);
        idNode->parent = this;
        this->list.push_back(std::move(idNode));
        state = State::NEXT;
      } else if (stream[0].type == TokenType::LBRACKET) {
        auto listNode = std::unique_ptr<ListNode>(new ListNode());
        listNode->parse(stream);
        listNode->parent = this;
        this->list.push_back(std::move(listNode));
        state = State::NEXT;
      } else if (stream[0].type == TokenType::LBRACE) {
        auto dictNode = std::unique_ptr<DictNode>(new DictNode());
        dictNode->parse(stream);
        dictNode->parent = this;
        this->list.push_back(std::move(dictNode));
        state = State::NEXT;
      } else {
        throw json_exception(
          stream[0], "expected list value type (string, id, list or dict)");
      }
    } else if (state == State::NEXT) {
      if (stream[0].type == TokenType::COMMA) {
        stream.erase(stream.begin());
        state = State::ITEMVALUE;
      } else if (stream[0].type == TokenType::RBRACKET) {
        stream.erase(stream.begin());
        return;
      } else {
        throw json_exception(stream[0], "expected \",\" or \"]\"");
      }
    }
  }

  throw json_exception("unexpected end-of-file");
}

Node& ListNode::operator[](size_t index) {
  return *this->list[index];
}

void ListNode::serialize(std::ofstream& outFile, size_t indentWidth) {
  outFile << "[";
  bool first = true;

  for (const std::unique_ptr<Node>& node : this->list) {
    if (first) {
      first = false;
    } else {
      outFile << ", ";
    }

    node->serialize(outFile, indentWidth + Node::SERIALIZE_INDENT);
  }

  outFile << "]";
}

size_t ListNode::size() {
  return this->list.size();
}

void ListNode::addValue(std::unique_ptr<Node> node) {
  if (node->parent != nullptr) {
    throw json_exception("addItem(): value was already added");
  }

  node->parent = this;
  this->list.push_back(std::move(node));
}

Node* ListNode::clone() {
  ListNode* newNode = new ListNode(*this);
  return newNode;
}

std::unique_ptr<Node> ListNode::removeValue(size_t index) {
  if (index >= list.size()) {
    throw json_exception("removeItem(): index is out-of-bounds");
  }

  auto node = std::move(this->list[index]);
  node->parent = nullptr;
  this->list.erase(this->list.begin() + index);
  return node;
}

// returns the list node to which the value was added
Node& ListNode::addTextValue(const std::string& value) {
  auto textNode = std::unique_ptr<TextNode>(new TextNode());
  textNode->set(value);
  this->addValue(std::move(textNode));
  return *this;
}

// returns the list node to which the value was added
Node& ListNode::addIdValue(const std::string& value) {
  auto idNode = std::unique_ptr<IDNode>(new IDNode());
  idNode->set(value);
  this->addValue(std::move(idNode));
  return *this;
}

// returns the list node to which the value was added
Node& ListNode::addIdValue(const char* value) {
  return this->addIdValue(std::string(value));
}

// returns the list node to which the value was added
Node& ListNode::addIdValue(const double& value) {
  auto idNode = std::unique_ptr<IDNode>(new IDNode());
  idNode->setDouble(value);
  this->addValue(std::move(idNode));
  return *this;
}

// returns the list node to which the value was added
Node& ListNode::addIdValue(const uint64_t& value) {
  auto idNode = std::unique_ptr<IDNode>(new IDNode());
  idNode->setUInt(value);
  this->addValue(std::move(idNode));
  return *this;
}

// returns the list node to which the value was added
Node& ListNode::addIdValue(const int64_t& value) {
  auto idNode = std::unique_ptr<IDNode>(new IDNode());
  idNode->setInt(value);
  this->addValue(std::move(idNode));
  return *this;
}

// returns the list node to which the value was added
Node& ListNode::addIdValue(const bool& value) {
  auto idNode = std::unique_ptr<IDNode>(new IDNode());
  idNode->setBool(value);
  this->addValue(std::move(idNode));
  return *this;
}

// returns created dict node
Node& ListNode::addDictValue() {
  auto dictNode = std::unique_ptr<DictNode>(new DictNode());
  auto& reference = *dictNode;  // because dictNode will be invalidated
  this->addValue(std::move(dictNode));
  return reference;
}

// returns created dict node
Node& ListNode::addListValue() {
  auto listNode = std::unique_ptr<ListNode>(new ListNode());
  auto& reference = *listNode;  // because listNode will be invalidated
  this->addValue(std::move(listNode));
  return reference;
}

std::unique_ptr<Node> ListNode::erase(Node& node) {
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

}  // namespace json
