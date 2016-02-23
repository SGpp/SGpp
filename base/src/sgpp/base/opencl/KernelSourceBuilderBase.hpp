// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "OCLDevice.hpp"
#include "sgpp/base/tools/json/Node.hpp"

namespace SGPP {
namespace base {

template <typename T>
class KernelSourceBuilderBase {
 protected:
  std::vector<std::string> indent;

  static const size_t MAX_INDENT_LEVEL = 10;

  //    std::string indent2;
  //
  //    std::string indent3;
  //
  //    std::string indent4;
  //
  //    std

  std::string floatType() {
    return std::is_same<T, float>::value ? "float" : "double";
  }

  std::string constSuffix() { return std::is_same<T, float>::value ? "f" : ""; }

  std::string intType() {
    return std::is_same<T, float>::value ? "uint" : "ulong";
  }

  std::string reuseSource(std::string fileName) {
    std::stringstream sourceStream;
    std::ifstream file;
    file.open(fileName);

    if (file.is_open()) {
      std::string line;

      while (getline(file, line)) {
        sourceStream << line << std::endl;
      }

      file.close();
    } else {
      throw base::operation_exception(
          "OCL error: file to reuse not found\n");
    }

    return sourceStream.str();
  }

  void writeSource(std::string fileName, std::string source) {
    // update file with kernel (for debugging)
    std::ofstream sourceFile;
    sourceFile.open(fileName);
    sourceFile << source;
    sourceFile.close();
  }

 public:
  KernelSourceBuilderBase() {
    //:
    //            indent(4, ' '), indent2(8, ' '), indent3(12, ' '), indent4(16,
    //            ' ') {
    for (size_t i = 1; i < MAX_INDENT_LEVEL + 1; i++) {
      indent.emplace_back(4 * i, ' ');
    }
  }
};
}  // namespace base
}  // namespace SGPP
