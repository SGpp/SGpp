// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef FUNCTIONLOOKUPTABLE_HPP_
#define FUNCTIONLOOKUPTABLE_HPP_

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/MultiFunction.hpp>

#include <vector>
#include <memory>
#include <unordered_map>
#include <string>
#include <mutex>

namespace sgpp {
namespace combigrid {

class DataVectorEqualTo {
 public:
  bool operator()(base::DataVector const &first, base::DataVector const &second) const {
    if (first.getSize() != second.getSize()) {
      return false;
    }

    for (size_t i = 0; i < first.getSize(); ++i) {
      if (first[i] != second[i]) {
        return false;
      }
    }

    return true;
  }
};

class DataVectorHash {
 public:
  size_t operator()(base::DataVector const &vec) const {
    std::hash<double> h;
    size_t result = 0;
    for (size_t i = 0; i < vec.getSize(); ++i) {
      result ^= h(vec[i]) * (i + 1);
    }
    return result;
  }
};

class FunctionLookupTable {
  std::shared_ptr<std::unordered_map<base::DataVector, double, DataVectorHash, DataVectorEqualTo>>
      hashmap;
  MultiFunction func;
  std::mutex tableMutex;

 public:
  FunctionLookupTable(MultiFunction func);

  double operator()(base::DataVector const &x);

  double eval(base::DataVector const &x);
  double evalThreadsafe(base::DataVector const &x);

  bool containsEntry(base::DataVector const &x);

  void addEntry(base::DataVector const &x, double y);

  std::string serialize();

  void deserialize(std::string const &value);

  size_t getNumEntries() const;
};
}  // namespace combigrid
}  // namespace sgpp

#endif /* FUNCTIONLOOKUPTABLE_HPP_ */
