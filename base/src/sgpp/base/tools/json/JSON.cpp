// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>

#include <fstream>
#include <sstream>
#include <vector>
#include <string>

namespace json {

JSON::JSON() : fileName("") {}

JSON::JSON(const std::string& fileName) : fileName(fileName) {
  std::ifstream file;
  file.exceptions(std::ifstream::failbit | std::ifstream::badbit);

  try {
    file.open(fileName);
  } catch (std::ifstream::failure&) {
    std::stringstream stream;
    stream << "json error: could not open file: " << fileName << std::endl;
    throw json_exception(stream.str());
  }

  std::string content;
  try {
    content.assign(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
  } catch (std::ifstream::failure&) {
    std::stringstream stream;
    stream << "json error: could not successfully read file: " << fileName << std::endl;
    throw json_exception(stream.str());
  }

  file.close();

  std::vector<Token> tokenStream = this->tokenize(content);

  this->parse(tokenStream);

  if (tokenStream.size() != 0) {
    throw json_exception(tokenStream[0], "expected end-of-file");
  }
}

JSON::JSON(const JSON& original) : DictNode(original) { this->fileName = original.fileName; }

JSON& JSON::operator=(const JSON& original) {
  DictNode::operator=(original);
  this->fileName = original.fileName; 
  return *this;
}

std::vector<Token> JSON::tokenize(const std::string& input) {
  std::vector<Token> stream;

  TokenType state = TokenType::NONE;

  Token token;
  size_t lineNumber = 1;
  size_t charNumber = 0;

  for (size_t i = 0; i < input.size(); i++) {
    if (input[i] == '\n') {
      lineNumber += 1;
      charNumber = 0;
    } else {
      charNumber += 1;
    }

    // skip whitespace while not tokenizing anything
    if (state == TokenType::NONE) {
      if (input[i] == ' ' || input[i] == '\r' || input[i] == '\n' || input[i] == '\t') {
        continue;
      }
    }

    if (state == TokenType::NONE) {
      if (input[i] == '{') {
        token.type = TokenType::LBRACE;
        token.value = "{";
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        stream.push_back(token);
      } else if (input[i] == ',') {
        token.type = TokenType::COMMA;
        token.value = ",";
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        stream.push_back(token);
      } else if (input[i] == '}') {
        token.type = TokenType::RBRACE;
        token.value = "}";
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        stream.push_back(token);
      } else if (input[i] == '[') {
        token.type = TokenType::LBRACKET;
        token.value = "[";
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        stream.push_back(token);
      } else if (input[i] == ']') {
        token.type = TokenType::RBRACKET;
        token.value = "]";
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        stream.push_back(token);
      } else if (input[i] == ':') {
        token.type = TokenType::COLON;
        token.value = ":";
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        stream.push_back(token);
      } else if (input[i] == '"') {
        token.type = TokenType::STRING;
        token.value = "";
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        state = TokenType::STRING;
      } else if (input[i] == '/') {
        state = TokenType::COMMENT;
      } else {
        token.type = TokenType::ID;
        token.value = input[i];
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        state = TokenType::ID;
      }
    } else if (state == TokenType::STRING) {
      if (input[i] == '"') {
        stream.push_back(token);
        state = TokenType::NONE;
      } else {
        token.value.push_back(input[i]);
      }
    } else if (state == TokenType::ID) {
      if (input[i] == '{' || input[i] == '}' || input[i] == '[' || input[i] == ']' ||
          input[i] == ':' || input[i] == ',' || input[i] == ' ' || input[i] == '\t' ||
          input[i] == '\r' || input[i] == '\n') {
        stream.push_back(token);
        state = TokenType::NONE;

        if (input[i] == '\n') {
          lineNumber -= 1;  // as the char will be reprocessed
        }

        i -= 1;  // revert by one

      } else {
        token.value.push_back(input[i]);
      }
    } else if (state == TokenType::COMMENT) {
      if (input[i] == '/') {
        state = TokenType::SINGLELINE;
      } else if (input[i] == '*') {
        state = TokenType::MULTILINECOMMENT;
      } else {
        std::stringstream messageStream;
        messageStream << "error: (line: " << lineNumber << ", char: " << (charNumber - 1)
                      << "): expected a single- or multiline comment after \"/\"";
        throw json_exception(messageStream.str());
      }
    } else if (state == TokenType::SINGLELINE) {
      if (input[i] == '\n') {
        state = TokenType::NONE;
      }
    } else if (state == TokenType::MULTILINECOMMENT) {
      if (input[i] == '*') {
        state = TokenType::MULTILINECOMMENTSTAR;
      }
    } else if (state == TokenType::MULTILINECOMMENTSTAR) {
      if (input[i] == '/') {
        state = TokenType::NONE;
      }
    } else {
      std::stringstream messageStream;
      messageStream << "error: (line: " << lineNumber << ", char: " << (charNumber - 1)
                    << "): illegal parser state, this might be a bug";
      throw json_exception(messageStream.str());
    }
  }

  return stream;
}

void JSON::serialize(const std::string& outFileName) {
  std::ofstream outFile(outFileName);
  outFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

  this->serialize(outFile, 0);

  outFile.close();
}

void JSON::deserialize(std::string content) {
  this->clear();
  std::vector<Token> tokenStream = this->tokenize(content);

  this->parse(tokenStream);

  if (tokenStream.size() != 0) {
    throw json_exception(tokenStream[0], "expected end-of-file");
  }
  fileName = "";
}
JSON* JSON::clone() { return new JSON(*this); }

void JSON::clear() {
  this->fileName = "";
  this->attributes.clear();
  this->keyOrder.clear();
  // it shouldn't even be necessary to reset these
  this->orderedKeyIndex = 0;
  this->parent = nullptr;
}

void JSON::deserializeFromString(const std::string &content) {
    this->clear();
    std::vector<Token> tokenStream = this->tokenize(content);

    this->parse(tokenStream);

    if (tokenStream.size() != 0) {
      throw json_exception(tokenStream[0], "expected end-of-file");
    }
}

}  // namespace json
