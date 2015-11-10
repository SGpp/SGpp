/*
 * ConfigurationParser.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#include <fstream>
#include <sstream>

#include "JSON.hpp"
#include "json_exception.hpp"

namespace json {

JSON::JSON(): fileName("") {

}

JSON::JSON(std::string fileName) :
    fileName(fileName) {
  std::ifstream file(fileName);
  file.exceptions(std::ifstream::failbit | std::ifstream::badbit);

  std::string content;
  content.assign(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());

  file.close();

  std::vector<JSONToken> tokenStream = this->tokenize(content);

  this->parse(tokenStream);
  if (tokenStream.size() != 0) {
    throw json_exception(tokenStream[0], "expected end-of-file");
  }

}

std::vector<JSONToken> JSON::tokenize(std::string &input) {
  std::vector<JSONToken> stream;

  JSONTokenType state = JSONTokenType::NONE;

  JSONToken token;
  size_t lineNumber = 1;
  size_t charNumber = 0;

  for (size_t i = 0; i < input.size(); i++) {
    if (input[i] == '\n') {
      lineNumber += 1;
      charNumber = 0;
    } else {
      charNumber += 1;
    }

//skip whitespace while not tokenizing anything
    if (state == JSONTokenType::NONE) {
      if (input[i] == ' ' || input[i] == '\r' || input[i] == '\n' || input[i] == '\t') {
        continue;
      }
    }

    if (state == JSONTokenType::NONE) {
      if (input[i] == '{') {
        token.type = JSONTokenType::LBRACE;
        token.value = "";
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        stream.push_back(token);
      } else if (input[i] == ',') {
        token.type = JSONTokenType::COMMA;
        token.value = "";
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        stream.push_back(token);
      } else if (input[i] == '}') {
        token.type = JSONTokenType::RBRACE;
        token.value = "";
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        stream.push_back(token);
      } else if (input[i] == '[') {
        token.type = JSONTokenType::LBRACKET;
        token.value = "";
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        stream.push_back(token);
      } else if (input[i] == ']') {
        token.type = JSONTokenType::RBRACKET;
        token.value = "";
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        stream.push_back(token);
      } else if (input[i] == ':') {
        token.type = JSONTokenType::COLON;
        token.value = "";
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        stream.push_back(token);
      } else if (input[i] == '"') {
        token.type = JSONTokenType::STRING;
        token.value = "";
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        state = JSONTokenType::STRING;
      } else if (input[i] == '/') {
        state = JSONTokenType::COMMENT;
      } else {
        token.type = JSONTokenType::ID;
        token.value = input[i];
        token.lineNumber = lineNumber;
        token.charNumber = charNumber;
        state = JSONTokenType::ID;
      }
    } else if (state == JSONTokenType::STRING) {
      if (input[i] == '"') {
        stream.push_back(token);
        state = JSONTokenType::NONE;
      } else {
        token.value.push_back(input[i]);
      }
    } else if (state == JSONTokenType::ID) {
      if (input[i] == '{' || input[i] == '}' || input[i] == '[' || input[i] == ']' || input[i] == ':' || input[i] == ','
          || input[i] == ' ' || input[i] == '\t' || input[i] == '\r' || input[i] == '\n') {
        stream.push_back(token);
        state = JSONTokenType::NONE;
        if (input[i] == '\n') {
          lineNumber -= 1; //as the char will be reprocessed
        }
        i -= 1; // revert by one

      } else {
        token.value.push_back(input[i]);
      }
    } else if (state == JSONTokenType::COMMENT) {
      if (input[i] == '/') {
        state = JSONTokenType::SINGLELINE;
      } else if (input[i] == '*') {
        state = JSONTokenType::MULTILINECOMMENT;
      } else {
        std::stringstream messageStream;
        messageStream << "error: (line: " << lineNumber << ", char: " << (charNumber - 1) << "): expected a single- or multiline comment after \"/\"";
        throw json_exception(messageStream.str());
      }
    } else if (state == JSONTokenType::SINGLELINE) {
      if (input[i] == '\n') {
        state = JSONTokenType::NONE;
      }
    } else if (state == JSONTokenType::MULTILINECOMMENT) {
      if (input[i] == '*') {
        state = JSONTokenType::MULTILINECOMMENTSTAR;
      }
    } else if (state == JSONTokenType::MULTILINECOMMENTSTAR) {
      if (input[i] == '/') {
        state = JSONTokenType::NONE;
      }
    } else {
      std::stringstream messageStream;
      messageStream << "error: (line: " << lineNumber << ", char: " << (charNumber - 1) << "): illegal parser state, this might be a bug";
      throw json_exception(messageStream.str());
    }
  }

  return stream;
}

void JSON::serialize(std::string outFileName) {
  std::ofstream outFile(outFileName);
  outFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

  this->serialize(outFile, 0);

  outFile.close();
}

}
