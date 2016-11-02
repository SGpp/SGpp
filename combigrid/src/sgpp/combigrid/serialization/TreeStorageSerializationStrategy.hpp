// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TREESTORAGESERIALIZATIONSTRATEGY_HPP_
#define TREESTORAGESERIALIZATIONSTRATEGY_HPP_

#include <sgpp/combigrid/common/MultiIndexIterator.hpp>
#include <sgpp/combigrid/serialization/AbstractSerializationStrategy.hpp>
#include <sgpp/combigrid/serialization/DefaultSerializationStrategy.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>
#include <sgpp/globaldef.hpp>

#include <stdexcept>
#include <string>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * This class provides a serialization strategy for TreeStorage<T>-objects (though not particularly
 * well-compressed).
 */
template <typename T>
class TreeStorageSerializationStrategy
    : public AbstractSerializationStrategy<std::shared_ptr<TreeStorage<T>>> {
  std::shared_ptr<AbstractSerializationStrategy<T>> innerStrategy;
  size_t numDimensions;

  static std::shared_ptr<AbstractSerializationStrategy<T>> getDefaultStrategy() {
    return std::shared_ptr<AbstractSerializationStrategy<T>>(new DefaultSerializationStrategy<T>());
  }

 public:
  /**
   * Constructor. A serialization strategy for the type T can be provided, otherwise, a
   * DefaultSerializationStrategy<T> will be used. The serialization strategy may use any
   * characters, this class will use escaping to deal with it.
   * @param numDimensions Dimension of the tree storage.
   */
  TreeStorageSerializationStrategy(
      size_t numDimensions,
      std::shared_ptr<AbstractSerializationStrategy<T>> innerStrategy = getDefaultStrategy())
      : innerStrategy(innerStrategy), numDimensions(numDimensions) {}

  virtual ~TreeStorageSerializationStrategy() {}

  virtual std::string serialize(std::shared_ptr<TreeStorage<T>> const &storage) {
    DefaultSerializationStrategy<size_t> indexStrategy;

    std::vector<std::string> entries;

    for (auto it = storage->getStoredDataIterator(); it->isValid(); it->moveToNext()) {
      std::vector<std::string> vectorEntries;

      auto multiIndex = it->getMultiIndex();

      for (size_t i = 0; i < multiIndex.size(); ++i) {
        vectorEntries.push_back(indexStrategy.serialize(multiIndex[i]));
      }

      entries.push_back(join(vectorEntries, ", ") + " -> " +
                        escape(innerStrategy->serialize(it->value()), '\\', ",->\n", "cmrn"));
    }

    return join(entries, "\n");
  }

  virtual std::shared_ptr<TreeStorage<T>> deserialize(std::string const &input) {
    std::shared_ptr<TreeStorage<T>> storage(new TreeStorage<T>(numDimensions));

    DefaultSerializationStrategy<size_t> indexStrategy;

    std::vector<std::string> entries = split(input, "\n");

    for (size_t i = 0; i < entries.size(); ++i) {
      std::vector<std::string> keyValuePair = split(entries[i], " -> ");

      if (keyValuePair.size() != 2) {
        throw std::runtime_error(
            "TreeStorageSerializationStrategy::deserialize(): keyValuePair does not have two "
            "entries");
      }

      std::vector<std::string> multiIndexEntries = split(keyValuePair[0], ", ");

      MultiIndex multiIndex(multiIndexEntries.size());

      for (size_t j = 0; j < multiIndexEntries.size(); ++j) {
        multiIndex[j] = indexStrategy.deserialize(multiIndexEntries[j]);
      }

      T value = innerStrategy->deserialize(unescape(keyValuePair[1], '\\', ",->\n", "cmrn"));

      storage->set(multiIndex, value);
    }

    return storage;
  }
};
}  // namespace combigrid
}  // namespace sgpp

#endif /* TREESTORAGESERIALIZATIONSTRATEGY_HPP_ */
