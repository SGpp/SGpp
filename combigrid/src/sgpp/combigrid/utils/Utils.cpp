// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/utils/Utils.hpp>

#include <string>
#include <vector>
#include <stdexcept>

namespace sgpp {
namespace combigrid {

std::int64_t binom(std::int64_t n, std::int64_t k) {
  if (k < 0 || n < k) {
    return 0;
  }

  if (2 * k > n) {
    k = n - k;
  }

  std::int64_t prod = 1;

  for (std::int64_t i = 1; i <= k; ++i) {
    prod *= (n + i - k);
    prod /= i;
  }

  return prod;
}

std::vector<std::string> split(std::string str, std::string separator) {
  std::vector<std::string> result;

  if (str.length() == 0) {
    return result;
  }

  size_t startPos = 0;

  while (true) {
    size_t afterEndPos = str.find(separator, startPos);

    if (afterEndPos == std::string::npos) {
      result.push_back(str.substr(startPos));
      return result;
    }

    result.push_back(str.substr(startPos, afterEndPos - startPos));
    startPos = afterEndPos + separator.length();
  }
}

std::string join(std::vector<std::string> const &elements, std::string const &separator) {
  if (elements.size() == 0) {
    return "";
  }

  std::string result = elements[0];

  for (size_t i = 1; i < elements.size(); ++i) {
    result += separator + elements[i];
  }

  return result;
}

std::string escape(std::string str, char escapeCharacter, std::string avoidCharacters,
                   std::string replaceCharacters) {
  if (avoidCharacters.length() != replaceCharacters.length()) {
    throw std::runtime_error("escape(): avoidCharacters.length() != replaceCharacters.length()");
  }

  std::string result;

  for (size_t i = 0; i < str.length(); ++i) {
    if (str[i] == escapeCharacter) {
      result += escapeCharacter;
      result += escapeCharacter;
    } else {
      size_t avoidPos = avoidCharacters.find(str[i]);

      if (avoidPos == std::string::npos) {
        result += str[i];
      } else {
        result += escapeCharacter;
        result += replaceCharacters[avoidPos];
      }
    }
  }

  return result;
}

std::string unescape(std::string str, char escapeCharacter, std::string avoidCharacters,
                     std::string replaceCharacters) {
  if (avoidCharacters.length() != replaceCharacters.length()) {
    throw std::runtime_error("unescape(): avoidCharacters.length() != replaceCharacters.length()");
  }

  std::string result;

  for (size_t i = 0; i < str.length(); ++i) {
    if (str[i] == escapeCharacter) {
      ++i;

      if (i >= str.length()) {
        throw std::runtime_error("unescape(): string ended after escape character");
      }

      char next = str[i];

      if (next == escapeCharacter) {
        result += escapeCharacter;
      } else {
        size_t replacePos = replaceCharacters.find(next);

        if (replacePos == std::string::npos) {
          throw std::runtime_error("unescape(): illegal escape character");
        } else {
          result += avoidCharacters[replacePos];
        }
      }
    } else {
      result += str[i];
    }
  }

  return result;
}
}  // namespace combigrid
}  // namespace sgpp
