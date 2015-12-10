/*
 * JSONDictNode.hpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#pragma once

#include <map>
#include <memory>

#include "Node.hpp"

namespace json {

  class DictNode: public Node {
    protected:
      std::map<std::string, std::unique_ptr<Node>> attributes;

      std::vector<std::string> keyOrder;

    public:

      DictNode();

      DictNode(const DictNode& original);

      DictNode& operator=(const DictNode& right);

      virtual void parse(std::vector<Token>& stream) override;

      void parseAttributes(std::vector<Token>& stream);

      virtual void serialize(std::ofstream& outFile, size_t indentWidth) override;

      virtual Node& operator[](const std::string& key) override;

      virtual size_t size() override;

      virtual Node* clone() override;

      virtual void addAttribute(const std::string& name, std::unique_ptr<Node> node) override;

      virtual std::unique_ptr<Node> removeAttribute(const std::string& name) override;

      // returns the node to which the attribute was added
      virtual Node& addTextAttr(const std::string& name, const std::string& value) override;

      // returns the node to which the attribute was added
      virtual Node& addIDAttr(const std::string& name, const std::string& value) override;

      // returns the node to which the attribute was added
      virtual Node& addIDAttr(const std::string& name, const double& value) override;

      // returns the node to which the attribute was added
      virtual Node& addIDAttr(const std::string& name, const uint64_t& value) override;

      // returns the node to which the attribute was added
      virtual Node& addIDAttr(const std::string& name, const int64_t& value) override;

      // returns the node to which the attribute was added
      virtual Node& addIDAttr(const std::string& name, const bool& value) override;

      // returns the node to which the attribute was added
      virtual Node& addIDAttr(const std::string& name, const char* value) override;

      // returns created dict node
      virtual Node& addDictAttr(const std::string& name) override;

      // returns created list node
      virtual Node& addListAttr(const std::string& name) override;

      // returns the node to which the attribute was added
      // replaces a node, adds a new node, if the node does not exist, the old node is deleted
      virtual Node& replaceTextAttr(const std::string& name, const std::string& value) override;

      // returns the node to which the attribute was added
      // replaces a node, adds a new node, if the node does not exist, the old node is deleted
      virtual Node& replaceIDAttr(const std::string& name, const std::string& value) override;

      // returns the node to which the attribute was added
      // replaces a node, adds a new node, if the node does not exist, the old node is deleted
      virtual Node& replaceIDAttr(const std::string& name, const double& value) override;

      // returns the node to which the attribute was added
      // replaces a node, adds a new node, if the node does not exist, the old node is deleted
      virtual Node& replaceIDAttr(const std::string& name, const uint64_t& value) override;

      // returns the node to which the attribute was added
      // replaces a node, adds a new node, if the node does not exist, the old node is deleted
      virtual Node& replaceIDAttr(const std::string& name, const int64_t& value) override;

      // returns the node to which the attribute was added
      // replaces a node, adds a new node, if the node does not exist, the old node is deleted
      virtual Node& replaceIDAttr(const std::string& name, const bool& value) override;

      // returns created dict node
      // replaces a node, adds a new node, if the node does not exist, the old node is deleted
      virtual Node& replaceDictAttr(const std::string& name) override;

      // returns created list node
      // replaces a node, adds a new node, if the node does not exist, the old node is deleted
      virtual Node& replaceListAttr(const std::string& name) override;

      virtual bool contains(const std::string& key) override;

      virtual std::unique_ptr<Node> erase(Node& node) override;

      virtual std::vector<std::string>& keys() override;
  };

}
