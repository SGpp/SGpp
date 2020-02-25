// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

class StringTokenizer {
 public:
  static void tokenize(const std::string &s, const char *delim,
                       std::vector<std::string> &v) {
    size_t pos = 0;
    size_t start = 0;

    // for each token
    while (start != std::string::npos) {
      // search for the first ocurence of the delimiter
      pos = s.find(delim, start);
      v.push_back(s.substr(start, (pos - start)));

      start = pos != std::string::npos ? ++pos : std::string::npos;
    }
  }
};
} /* namespace base */
} /* namespace sgpp */
