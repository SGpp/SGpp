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

  // this checks if level (1,...,1) exists. due to the different point generation on level 0 and 1
  // in the combigrid module and classical SG++ the combigrid grid has no exactly matching SG++ grid
  // (?)
  //  bool containsLevelOne = false;
  //  MultiIndex levelOne(d, 1);
  //  while (it->isValid()) {
  //    MultiIndex currentLevel = it->getMultiIndex();
  //    if (currentLevel == levelOne) {
  //      containsLevelOne = true;
  //    }
  //    it->moveToNext();
  //  }
  //  if (!containsLevelOne) {
  //    std::cerr
  //        << "GridConversion Warning: The subgrid of level 1 in every dimension (1,...,1) must
  //        exist."
  //        << std::endl;
  //  }

  it = levelStructure->getStoredDataIterator();
  while (it->isValid()) {
    MultiIndex currentLevel = it->getMultiIndex();

    // now iterate over index
    MultiIndex multiBounds(d);
    for (size_t i = 0; i < d; ++i) {
      // there are 2^l+1 elements in the level
      multiBounds[i] = (1 << currentLevel[i]) + 1;
    }
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
          isValidPoint = false;
          break;
        }
        point.push(i, sglevel, sgindex);
      }
      if (isValidPoint) {
        point.rehash();
        if (!storage.isContaining(point)) {
          storage.insert(point);
        }
      }
      indexIt.moveToNext();
    }

    it->moveToNext();
  }

  storage.recalcLeafProperty();
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

sgpp::base::DataMatrix convertLevelStructureToGridPoints(
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const& levelStructure,
    size_t numDimensions) {
  sgpp::base::DataMatrix gridpointMatrix(0, numDimensions);
  sgpp::base::GridStorage gridStorage(numDimensions);
  convertexpUniformBoundaryCombigridToHierarchicalSparseGrid(levelStructure, gridStorage);
  for (size_t q = 0; q < gridStorage.getSize(); q++) {
    auto point = gridStorage.getPoint(q);
    sgpp::base::DataVector row;
    for (size_t d = 0; d < gridStorage.getDimension(); d++) {
      row.push_back(point.getStandardCoordinate(d));
    }
    gridpointMatrix.appendRow(row);
  }
  return gridpointMatrix;
}

sgpp::base::DataVector calculateInterpolationCoefficientsForConvertedExpUniformBoundaryCombigird(
    std::shared_ptr<sgpp::base::Grid>& grid, sgpp::base::GridStorage& gridStorage,
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation>& combigridInterpolationOperation,
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const& levelStructure) {
  sgpp::base::DataMatrix interpolParams(gridStorage.getDimension(), gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
      p[j] = gp.getStandardCoordinate(j);
    }
    interpolParams.setColumn(i, p);
  }

  // obtain function values from combigrid surrogate
  combigridInterpolationOperation->setParameters(interpolParams);

  combigridInterpolationOperation->getLevelManager()->addLevelsFromStructure(levelStructure);
  sgpp::base::DataVector f_values = combigridInterpolationOperation->getResult();

  sgpp::base::Printer::getInstance().setVerbosity(-1);
  sgpp::base::HierarchisationSLE hierSLE(*grid);
  sgpp::base::sle_solver::Auto sleSolver;
  sgpp::base::DataVector alpha(grid->getSize());
  if (!sleSolver.solve(hierSLE, f_values, alpha)) {
    std::cout << "Solving failed!" << std::endl;
  }
  return alpha;
}

} /* namespace combigrid */
} /* namespace sgpp */
