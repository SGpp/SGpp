// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/combigrid/common/GridConversion.hpp>

#include <numeric>

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

// a SG grid of level 0 contains all vertices of the d dim hypercube. In combigrids these vertices
// appear the first time on level (1,...,1). Therefore this level must be contained in the combigrid
// to convert it to an SG grid.
void convertexpUniformBoundaryCombigridToHierarchicalSparseGrid(
    std::shared_ptr<TreeStorage<uint8_t>> levelStructure, base::HashGridStorage& storage) {
  size_t d = levelStructure->getNumDimensions();

  auto it = levelStructure->getStoredDataIterator();
  base::HashGridPoint::level_type sglevel(0);
  base::HashGridPoint::index_type sgindex(0);
  base::HashGridPoint point(d);

  // check if all level 1 grids exist, otherwise the grid cannot be transformed (at least not
  // without much more special cases)
  bool containsLevelOne = false;
  while (it->isValid()) {
    MultiIndex currentLevel = it->getMultiIndex();
    size_t levelsum = std::accumulate(currentLevel.begin(), currentLevel.end(), 0);
    if (levelsum == d) {
      containsLevelOne = true;
    }
    it->moveToNext();
  }
  if (!containsLevelOne) {
    std::cerr << "GridConversion: The subgrid of level 1 in every dimension (1,...,1) must exist."
              << std::endl;
  }

  it = levelStructure->getStoredDataIterator();
  while (it->isValid()) {
    MultiIndex currentLevel = it->getMultiIndex();

    //    std::cout << it->getMultiIndex()[0] << " ";

    // now iterate over index
    MultiIndex multiBounds(d);
    for (size_t i = 0; i < d; ++i) {
      // there are 2^l+1 elements in the level
      multiBounds[i] = (1 << currentLevel[i]) + 1;
      //      std::cout << multiBounds[i] << " ";
    }
    //    std::cout << "\n";
    MultiIndexIterator indexIt(multiBounds);
    while (indexIt.isValid()) {
      MultiIndex index = indexIt.getMultiIndex();
      bool isValidPoint = true;
      for (size_t i = 0; i < d; ++i) {
        if (currentLevel[i] == 0) {
          // point 0.5
          sglevel = static_cast<base::HashGridPoint::level_type>(1);
          sgindex = static_cast<base::HashGridPoint::index_type>(1);
        } else if (currentLevel[i] == 1) {
          // points 0, 0.5, 1
          if (index[i] == 0) {
            sglevel = static_cast<base::HashGridPoint::level_type>(0);
            sgindex = static_cast<base::HashGridPoint::index_type>(0);
          } else if (index[i] == 1) {
            sglevel = static_cast<base::HashGridPoint::level_type>(1);
            sgindex = static_cast<base::HashGridPoint::index_type>(1);
          } else if (index[i] == 2) {
            sglevel = static_cast<base::HashGridPoint::level_type>(0);
            sgindex = static_cast<base::HashGridPoint::index_type>(1);
          }
        } else if (index[i] % 2 == 1) {
          sglevel = static_cast<base::HashGridPoint::level_type>(currentLevel[i]);
          sgindex = static_cast<base::HashGridPoint::index_type>(index[i]);
        } else {
          //          std::cout << "X " << currentLevel[0] << " " << currentLevel[1] << " " <<
          //          index[0] << " "
          //                    << index[1] << std::endl;
          isValidPoint = false;
          break;
        }
        point.push(i, sglevel, sgindex);
      }
      if (isValidPoint) {
        point.rehash();
        if (!storage.isContaining(point)) {
          //          std::cout << point.getLevel(0) << " " << point.getLevel(1) << " " <<
          //          point.getIndex(0)
          //                    << " " << point.getIndex(1) << std::endl;
          storage.insert(point);
        }
      }
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
