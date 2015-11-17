/*
 * JSONIDNode.cpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#include "IDNode.hpp"

#include <fstream>
#include <string>
#include <sstream>

#include "json_exception.hpp"

namespace json {

IDNode::IDNode() :
        value(), internalType(InternalIDType::ID), doubleValue(0.0), unsignedValue(0), signedValue(0), boolValue(
                false) {
}

void IDNode::parse(std::vector<Token> &stream) {
//create new text node
    if (stream[0].type == TokenType::ID) {
        this->value = stream[0].value;
        stream.erase(stream.begin());

        this->setupInternalType();
    } else {
        throw json_exception(stream[0], "expected id");
    }
}

void IDNode::setupInternalType() {
    this->internalType = InternalIDType::ID;

    //try validating as bool
    if (this->value.compare("true") == 0) {
        this->internalType = InternalIDType::BOOL;
        this->boolValue = true;
        return;
    } else if (this->value.compare("false") == 0) {
        this->internalType = InternalIDType::BOOL;
        this->boolValue = false;
        return;
    }

    //try validating as unsigned integer
    try {
        std::string::size_type size;
        uint64_t asUnsigned = stoull(this->value, &size);
        if (this->value.size() == size) {
            this->unsignedValue = asUnsigned;
            this->internalType = InternalIDType::UINT;
            return;
        }
    } catch (std::invalid_argument &e) {
    }

    //try validating as signed integer
    try {
        std::string::size_type size;
        int64_t asSigned = stoll(this->value, &size);
        if (this->value.size() == size) {
            this->signedValue = asSigned;
            this->internalType = InternalIDType::INT;
            return;
        }
    } catch (std::invalid_argument &e) {
    }

    //try validating as double
    try {
        std::string::size_type size;
        double asDouble = stod(this->value, &size);
        if (this->value.size() == size) {
            this->doubleValue = asDouble;
            this->internalType = InternalIDType::DOUBLE;
            return;
        }
    } catch (std::invalid_argument &e) {
    }
}

std::string &IDNode::get() {
    return this->value;
}

void IDNode::set(const std::string &value) {
    this->value = value;

    this->setupInternalType();
}

double IDNode::getDouble() {
    if (this->internalType == InternalIDType::DOUBLE) {
        return this->doubleValue;
    } else {
        throw json_exception("node has not a numerical value");
    }
}

void IDNode::setDouble(double numericValue) {
    this->doubleValue = numericValue;
    this->internalType = InternalIDType::DOUBLE;
    std::stringstream stringstream;
    stringstream << numericValue;
    this->value = stringstream.str();
}

uint64_t IDNode::getUInt() {
    if (this->internalType == InternalIDType::UINT) {
        return this->unsignedValue;
    } else {
        throw json_exception("node has not a unsigned integer value");
    }
}

void IDNode::setUInt(uint64_t uintValue) {
    this->unsignedValue = uintValue;
    this->internalType = InternalIDType::UINT;
    std::stringstream stringstream;
    stringstream << uintValue;
    this->value = stringstream.str();
}

int64_t IDNode::getInt() {
    if (this->internalType == InternalIDType::INT) {
        return this->signedValue;
    } else {
        throw json_exception("node has not a unsigned integer value");
    }
}

void IDNode::setInt(int64_t intValue) {
    this->unsignedValue = intValue;
    this->internalType = InternalIDType::INT;
    std::stringstream stringstream;
    stringstream << intValue;
    this->value = stringstream.str();
}

bool IDNode::getBool() {
    if (this->internalType == InternalIDType::BOOL) {
        return this->boolValue;
    } else {
        throw json_exception("node has not a unsigned integer value");
    }
}

void IDNode::setBool(bool boolValue) {
    this->boolValue = boolValue;
    this->internalType = InternalIDType::BOOL;
    std::stringstream stringstream;
    stringstream << boolValue;
    this->value = stringstream.str();
}

void IDNode::serialize(std::ofstream &outFile, size_t indentWidth) {
    outFile << this->value;
}

size_t IDNode::size() {
    return 1;
}

Node *IDNode::clone() {
    IDNode *newNode = new IDNode(*this);
    return newNode;
}

}
