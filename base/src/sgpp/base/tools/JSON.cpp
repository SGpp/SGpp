/*
 * ConfigurationParser.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#include <fstream>
#include "JSON.hpp"

namespace SGPP {
namespace base {

JSON::JSON(std::string fileName) {
  std::ifstream file(fileName);

  if (!file.is_open()) {
    throw;
  }

//  std::string content(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());

  std::string content;
  content.assign(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());

//  std::string content("");

  file.close();

  std::vector<JSONToken> tokenStream = this->tokenize(content);

//  for (JSONToken &token : tokenStream) {
//    std::cout << "type: ";
//    if (token.type == JSONTokenType::COLON) {
//      std::cout << "COLON";
//    } else if (token.type == JSONTokenType::ID) {
//      std::cout << "ID value: " << token.value;
//    } else if (token.type == JSONTokenType::LBRACE) {
//      std::cout << "LBRACE";
//    } else if (token.type == JSONTokenType::LBRACKET) {
//      std::cout << "LBRACKET";
//    } else if (token.type == JSONTokenType::RBRACE) {
//      std::cout << "RBRACE";
//    } else if (token.type == JSONTokenType::RBRACKET) {
//      std::cout << "RBRACE";
//    } else if (token.type == JSONTokenType::STRING) {
//      std::cout << "STRING value: " << token.value;
//    } else if (token.type == JSONTokenType::COMMA) {
//      std::cout << "COMMA";
//    }
//    std::cout << std::endl;
//  }

  this->parse(tokenStream);
  if (tokenStream.size() != 0) {
    throw;
  }

}

std::vector<JSONToken> JSON::tokenize(std::string &input) {
  std::vector<JSONToken> stream;

  JSONTokenType state = JSONTokenType::NONE;

  JSONToken token;

  for (size_t i = 0; i < input.size(); i++) {
    //skip whitespace
    if (state != JSONTokenType::STRING) {
      if (input[i] == ' ' || input[i] == '\r' || input[i] == '\n' || input[i] == '\t') {
        continue;
      }
    }

    if (state == JSONTokenType::NONE) {
      if (input[i] == '{') {
        token.type = JSONTokenType::LBRACE;
        token.value = "";
        stream.push_back(token);
      } else if (input[i] == ',') {
        token.type = JSONTokenType::COMMA;
        token.value = "";
        stream.push_back(token);
      } else if (input[i] == '}') {
        token.type = JSONTokenType::RBRACE;
        token.value = "";
        stream.push_back(token);
      } else if (input[i] == '[') {
        token.type = JSONTokenType::LBRACKET;
        token.value = "";
        stream.push_back(token);
      } else if (input[i] == ']') {
        token.type = JSONTokenType::RBRACKET;
        token.value = "";
        stream.push_back(token);
      } else if (input[i] == ':') {
        token.type = JSONTokenType::COLON;
        token.value = "";
        stream.push_back(token);
      } else if (input[i] == '"') {
        token.type = JSONTokenType::STRING;
        token.value = "";
        state = JSONTokenType::STRING;
      } else {
        token.type = JSONTokenType::ID;
        token.value = input[i];
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
      if (input[i] == '{' || input[i] == '}' || input[i] == '[' || input[i] == ']' || input[i] == ':'
          || input[i] == ',') {
        stream.push_back(token);
        state = JSONTokenType::NONE;
        i -= 1; // revert by one
      } else {
        token.value.push_back(input[i]);
      }
    }
  }

  return stream;
}

void JSON::serialize(std::string outFileName) {
  std::ofstream outFile(outFileName);

  if (!outFile.is_open()) {
    throw;
  }

  this->serialize(outFile, 0);

  outFile.close();
}

//;
//}

}
}

