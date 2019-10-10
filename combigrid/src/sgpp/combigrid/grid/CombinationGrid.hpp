// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/basis/HeterogeneousBasis.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class CombinationGrid {
 public:
  CombinationGrid();

  CombinationGrid(const std::vector<FullGrid>& fullGrids, const base::DataVector& coefficients);

  static std::vector<LevelVector> enumerateLevelsWithSumWithBoundary(size_t dim, level_t n);

  static std::vector<LevelVector> enumerateLevelsWithSumWithoutBoundary(size_t dim, level_t n);

  static CombinationGrid fromRegular(size_t dim, level_t n, const HeterogeneousBasis& basis,
      bool hasBoundary = true);

  static CombinationGrid fromSubspaces(const std::vector<LevelVector>& subspaceLevels,
      const HeterogeneousBasis& basis, bool hasBoundary = true);

  void combinePoints(base::GridStorage& gridStorage) const;

  void combineSparseGridValues(const base::GridStorage& gridStorage,
      const std::vector<base::DataVector>& values, base::DataVector& result) const;

  void combineSparseGridValues(const base::GridStorage& gridStorage,
      const std::vector<base::DataMatrix>& values, base::DataMatrix& result) const;

  double combineValues(const base::DataVector& values) const;

  void combineValues(const base::DataMatrix& values, base::DataVector& result) const;

  void distributeValuesToFullGrid(const base::GridStorage& gridStorage,
      const base::DataVector& values, const FullGrid& fullGrid, base::DataVector& result) const;

  void distributeValuesToFullGrids(const base::GridStorage& gridStorage,
      const base::DataVector& values, std::vector<base::DataVector>& result) const;

  size_t getDimension() const;

  const std::vector<FullGrid>& getFullGrids() const;

  void setFullGrids(const std::vector<FullGrid>& fullGrids);

  const base::DataVector& getCoefficients() const;

  void setCoefficients(const base::DataVector& coefficients);

 protected:
  std::vector<FullGrid> fullGrids;
  base::DataVector coefficients;

  bool findGridPointInFullGrid(const FullGrid& fullGrid, const base::GridPoint& gridPoint,
      IndexVector& index) const;
};

}  // namespace combigrid
}  // namespace sgpp
