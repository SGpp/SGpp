// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <string>
#include <vector>

/**
 * This header contains some utility functions.
 */

namespace sgpp {
namespace combigrid {

/**
 * Exponentiation function for integer types with exact precision. Uses square-and-multiply.
 */
template <typename T>
T pow(T base, size_t exponent) {
  T result = 1;

  size_t numBits = 8 * sizeof(size_t);
  size_t mask = static_cast<size_t>(1) << (numBits - 1);

  // while (exponent != 0) {
  for (size_t i = 0; i < numBits; ++i) {
    result *= result;

    if (exponent & mask) {
      result *= base;
    }

    exponent = exponent << 1;
  }

  return result;
}

/**
 * Returns the binomial coefficient n over k.
 */
std::int64_t binom(std::int64_t n, std::int64_t k);

/**
 * Splits the string str at every occurrence of separator and returns the parts (without the
 * separator) in a vector.
 */
std::vector<std::string> split(std::string str, std::string separator);

/**
 * Concatenates the strings in elements and inserts the separator at the concatenation points. This
 * operation is inverse to split when the same separator is used.
 */
std::string join(std::vector<std::string> const &elements, std::string const &separator);

/**
 * Escapes in str each occurrence of a character in avoidCharacters with the escape character and
 * the corresponding character in replaceCharacters.
 * avoidCharacters and replaceCharacters must have the same length and use disjoint sets of
 * characters. The escapeCharacter is replaced with twice itself.
 * The escape character must differ from all characters in avoidCharacters and in replaceCharacters.
 * The resulting string will not contain any of the characters in avoidCharacters, which can be
 * useful for serialization.
 */
std::string escape(std::string str, char escapeCharacter, std::string avoidCharacters,
                   std::string replaceCharacters);

/**
 * Reverses the effect of the function escape(). The parameters have to fulfill the same
 * requirements as in escape() and the string str should have been generated with escape().
 */
std::string unescape(std::string str, char escapeCharacter, std::string avoidCharacters,
                     std::string replaceCharacters);

/**
 * Reads a file into a string (without advanced error-handling).
 */
std::string readFromFile(std::string filename);

/**
 * Writes a string into a file, overwriting currently saved data if the file already exists (without
 * advanced error-handling).
 */
void writeToFile(std::string filename, std::string value);
}  // namespace combigrid
}  // namespace sgpp

#endif /* UTILS_HPP_ */
