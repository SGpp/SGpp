/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * StringTokenizer.hpp
 *
 *  Created on: 21.06.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <cstring>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

// TODO (lettrich): Make the string tokenizer a bit more fancy.
class StringTokenizer {
 public:
  static void tokenize(const std::string& s, const char* delim, std::vector<std::string>& v) {
    // to avoid modifying original string
    // first duplicate the original string and return a char pointer then free the memory
    char* dup = strdup(s.c_str());
    char* token = strtok(dup, delim);
    while (token != NULL) {
      v.push_back(std::string(token));
      // the call is treated as a subsequent calls to strtok:
      // the function continues from where it left in previous invocation
      token = strtok(NULL, delim);
    }
    free(dup);
  };
};
} /* namespace datadriven */
} /* namespace sgpp */
