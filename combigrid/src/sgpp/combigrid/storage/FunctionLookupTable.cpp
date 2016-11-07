// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/serialization/DefaultSerializationStrategy.hpp>
#include <sgpp/combigrid/serialization/FloatSerializationStrategy.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

#include <stdexcept>
#include <string>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Helper class used internally as an equality predicate.
 */
class DataVectorEqualTo {
 public:
  bool operator()(base::DataVector const& first, base::DataVector const& second) const {
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

/**
 * Helper class used internally as a hash function for DataVector objects.
 */
class DataVectorHash {
 public:
  size_t operator()(base::DataVector const& vec) const {
    std::hash<double> h;
    size_t result = 0;
    for (size_t i = 0; i < vec.getSize(); ++i) {
      result ^= h(vec[i]) * (i + 1);
    }
    return result;
  }
};

/**
 * Helper to realize the PIMPL pattern
 */
struct FunctionLookupTableImpl {
  std::shared_ptr<std::unordered_map<base::DataVector, double, DataVectorHash, DataVectorEqualTo>>
      hashmap;
  MultiFunction func;
  std::mutex tableMutex;

  explicit FunctionLookupTableImpl(MultiFunction func)
      : hashmap(
            new std::unordered_map<base::DataVector, double, DataVectorHash, DataVectorEqualTo>()),
        func(func),
        tableMutex() {}
};

FunctionLookupTable::FunctionLookupTable(MultiFunction const& func)
    : impl(std::make_shared<FunctionLookupTableImpl>(func)) {}

double FunctionLookupTable::operator()(const base::DataVector& x) {
  auto it = impl->hashmap->find(x);
  if (it == impl->hashmap->end()) {
    auto y = impl->func(x);
    addEntry(x, y);
    return y;
  }

  return it->second;
}

double FunctionLookupTable::eval(const base::DataVector& x) { return (*this)(x); }

double FunctionLookupTable::evalThreadsafe(const base::DataVector& x) {
  impl->tableMutex.lock();
  auto it = impl->hashmap->find(x);
  if (it == impl->hashmap->end()) {
    impl->tableMutex.unlock();
    auto y = impl->func(x);
    impl->tableMutex.lock();
    addEntry(x, y);
    impl->tableMutex.unlock();
    return y;
  }

  std::lock_guard<std::mutex> guard(impl->tableMutex);
  auto y = it->second;
  impl->tableMutex.unlock();
  return y;
}

void FunctionLookupTable::addEntry(const base::DataVector& x, double y) { (*impl->hashmap)[x] = y; }

std::string FunctionLookupTable::serialize() {
  FloatSerializationStrategy<double> strategy;

  std::vector<std::string> entries;

  for (auto it = impl->hashmap->begin(); it != impl->hashmap->end(); ++it) {
    std::vector<std::string> vectorEntries;

    auto& vec = it->first;

    for (size_t i = 0; i < vec.getSize(); ++i) {
      vectorEntries.push_back(strategy.serialize(vec[i]));
    }

    entries.push_back(join(vectorEntries, ", ") + " -> " + strategy.serialize(it->second));
  }

  return join(entries, "\n");
}

void FunctionLookupTable::deserialize(const std::string& value) {
  FloatSerializationStrategy<double> strategy;

  std::vector<std::string> entries = split(value, "\n");

  for (size_t i = 0; i < entries.size(); ++i) {
    if (entries[i].length() == 0) {
      continue;
    }

    std::vector<std::string> keyValuePair = split(entries[i], " -> ");

    if (keyValuePair.size() != 2) {
      throw std::runtime_error(
          "FunctionLookupTable::deserialize(): keyValuePair does not have two entries but " +
          DefaultSerializationStrategy<size_t>().serialize(keyValuePair.size()) +
          " entries in line " + DefaultSerializationStrategy<size_t>().serialize(i) + ": " +
          entries[i]);
    }

    std::vector<std::string> dataVectorEntries = split(keyValuePair[0], ", ");

    base::DataVector x(dataVectorEntries.size());

    for (size_t j = 0; j < dataVectorEntries.size(); ++j) {
      x[j] = strategy.deserialize(dataVectorEntries[j]);
    }

    double y = strategy.deserialize(keyValuePair[1]);

    addEntry(x, y);
  }
}

bool FunctionLookupTable::containsEntry(const base::DataVector& x) {
  auto it = impl->hashmap->find(x);
  return it != impl->hashmap->end();
}

size_t FunctionLookupTable::getNumEntries() const { return impl->hashmap->size(); }

MultiFunction FunctionLookupTable::toMultiFunction() const { return MultiFunction(*this); }

}  // namespace combigrid
}  // namespace sgpp
