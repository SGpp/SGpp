// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/common/GridConversion.hpp>

namespace sgpp {
namespace combigrid {

std::shared_ptr<TreeStorage<uint8_t>> allStorageLevels(base::HashGridStorage* storage) {
  size_t d = storage->getDimension();
  auto treeStorage = std::make_shared<TreeStorage<uint8_t>>(d);

  auto it = storage->begin();

  for (; it != storage->end(); ++it) {
    MultiIndex level(d);
    for (size_t i = 0; i < d; ++i) {
      level[i] = it->first->getLevel(i) - 1;
    }

    treeStorage->set(level, 1);
  }

  return treeStorage;
}

base::HashGridStorage* toHashGridStorage(std::shared_ptr<TreeStorage<uint8_t>> levelStructure) {
  size_t d = levelStructure->getNumDimensions();
  auto storage = new base::HashGridStorage(d);

  auto it = levelStructure->getStoredDataIterator();

  while (it->isValid()) {
    MultiIndex currentLevel = it->getMultiIndex();

    // now iterate over index
    MultiIndex multiBounds(d);
    for (size_t i = 0; i < d; ++i) {
      // there are 2^l elements in the level, subtract 1 because the bound is inclusive
      multiBounds[i] = (1 << currentLevel[i]) - 1;
    }
    MultiIndexIterator indexIt(multiBounds);

    while (indexIt.isValid()) {
      MultiIndex index = indexIt.getMultiIndex();
      base::HashGridPoint point(d);

      for (size_t i = 0; i < d; ++i) {
        point.push(i, static_cast<base::HashGridPoint::level_type>(currentLevel[d] + 1),
                   static_cast<base::HashGridPoint::level_type>(2 * index[d] + 1));
      }

      point.rehash();

      storage->insert(point);

      indexIt.moveToNext();
    }

    it->moveToNext();
  }

  return storage;
}

} /* namespace combigrid */
} /* namespace sgpp */
