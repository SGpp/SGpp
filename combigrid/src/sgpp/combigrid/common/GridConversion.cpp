// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/combigrid/common/GridConversion.hpp>

namespace sgpp {
namespace combigrid {

std::shared_ptr<TreeStorage<uint8_t>> convertHierarchicalSparseGridToCombigrid(
    base::HashGridStorage& storage, GridConversionTypes conversionType) {
  switch (conversionType) {
    case GridConversionTypes::ALLSUBSPACES:
      return allStorageLevels(storage);
    case GridConversionTypes::COMPLETESUBSPACES:
    default:
      return completeStorageLevels(storage);
  }
}

void convertCombigridToHierarchicalSparseGrid(std::shared_ptr<TreeStorage<uint8_t>> levelStructure,
                                              base::HashGridStorage& storage) {
  size_t d = levelStructure->getNumDimensions();

  auto it = levelStructure->getStoredDataIterator();
  base::HashGridPoint::level_type sglevel;
  base::HashGridPoint::index_type sgindex;
  base::HashGridPoint point(d);

  while (it->isValid()) {
    MultiIndex currentLevel = it->getMultiIndex();
    std::cout << it->getMultiIndex()[0] << " ";

    // now iterate over index
    MultiIndex multiBounds(d);
    for (size_t i = 0; i < d; ++i) {
      // there are 2^l elements in the level, subtract 1 because the bound is inclusive
      multiBounds[i] = 1 << currentLevel[i];
      std::cout << multiBounds[i] << " ";
    }
    std::cout << "\n";
    MultiIndexIterator indexIt(multiBounds);

    while (indexIt.isValid()) {
      MultiIndex index = indexIt.getMultiIndex();
      for (size_t i = 0; i < d; ++i) {
        sglevel = static_cast<base::HashGridPoint::level_type>(currentLevel[i] + 1);
        sgindex = static_cast<base::HashGridPoint::index_type>(2 * index[i] + 1);
        point.push(i, sglevel, sgindex);
      }

      storage.insert(point);
      indexIt.moveToNext();
    }

    it->moveToNext();
  }

  storage.recalcLeafProperty();
}

// in combigrid a boundary grid of level 0 contains 0 .5,  a SG++ boundary grid of level 0 contains
// 0 and 1. From level 1 onwards they are identical. Therefore a combigrid boudnary grid of level 0
// cannot be converted into a SG++ boundary grid!
void convertBoundaryCombigridToHierarchicalSparseGrid(
    std::shared_ptr<TreeStorage<uint8_t>> levelStructure, base::HashGridStorage& storage) {
  size_t d = levelStructure->getNumDimensions();

  auto it = levelStructure->getStoredDataIterator();
  base::HashGridPoint::level_type sglevel;
  base::HashGridPoint::index_type sgindex;
  base::HashGridPoint point(d);

  // check if all level 1 grids exist, otherwise the grid cannot be transformed (at least not
  // without much more special cases)
  size_t numLevelOneLevels = 0;
  while (it->isValid()) {
    MultiIndex currentLevel = it->getMultiIndex();
    size_t levelsum = std::accumulate(currentLevel.begin(), currentLevel.end(), 0);
    if (levelsum == 1) {
      numLevelOneLevels++;
    }
  }
  if (numLevelOneLevels != d) {
    std::cerr << "GridConversion: The subgrid of level 1 must exist in every dimension."
              << std::endl;
  }

  it.reset();
  while (it->isValid()) {
    MultiIndex currentLevel = it->getMultiIndex();

    std::cout << it->getMultiIndex()[0] << " ";

    // now iterate over index
    MultiIndex multiBounds(d);
    for (size_t i = 0; i < d; ++i) {
      // there are 2^l elements in the level, subtract 1 because the bound is inclusive
      multiBounds[i] = 1 << currentLevel[i];
      std::cout << multiBounds[i] << " ";
    }
    std::cout << "\n";
    MultiIndexIterator indexIt(multiBounds);

    while (indexIt.isValid()) {
      MultiIndex index = indexIt.getMultiIndex();
      for (size_t i = 0; i < d; ++i) {
        // if level == 0 oder 1 Spezial Level und Index
        sglevel = static_cast<base::HashGridPoint::level_type>(currentLevel[i]);
        sgindex = static_cast<base::HashGridPoint::index_type>(2 * index[i] + 1);
        point.push(i, sglevel, sgindex);
      }

      storage.insert(point);
      indexIt.moveToNext();
    }

    it->moveToNext();
  }

  storage.recalcLeafProperty();

  // Funktionswerte berechnen via Interpolant aus Operation, Punkte aus HashGridStorage
}

std::shared_ptr<TreeStorage<uint8_t>> allStorageLevels(base::HashGridStorage& storage) {
  size_t numDims = storage.getDimension();
  auto treeStorage = std::make_shared<TreeStorage<uint8_t>>(numDims);

  auto it = storage.begin();

  for (; it != storage.end(); ++it) {
    MultiIndex level(numDims);
    for (size_t i = 0; i < numDims; ++i) {
      level[i] = it->first->getLevel(i) - 1;
    }

    treeStorage->set(level, 1);
  }

  return treeStorage;
}

std::shared_ptr<TreeStorage<uint8_t>> completeStorageLevels(base::HashGridStorage& storage) {
  size_t numDims = storage.getDimension();
  auto treeStorage = std::make_shared<TreeStorage<uint8_t>>(numDims);
  auto countStorage = std::make_shared<TreeStorage<size_t>>(numDims);

  auto it = storage.begin();

  for (; it != storage.end(); ++it) {
    MultiIndex level(numDims);
    for (size_t i = 0; i < numDims; ++i) {
      level[i] = it->first->getLevel(i) - 1;
    }

    ++countStorage->get(level);
  }

  auto countIt = countStorage->getStoredDataIterator();

  while (countIt->isValid()) {
    auto level = countIt->getMultiIndex();

    // compute number of points that are in the level if it is full
    size_t prod = 1;
    for (size_t dim = 0; dim < numDims; ++dim) {
      prod *= 1 << level[dim];
    }

    if (countIt->value() == prod) {
      treeStorage->set(level, 1);
    }
    countIt->moveToNext();
  }

  return treeStorage;
}

} /* namespace combigrid */
} /* namespace sgpp */
