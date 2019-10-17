// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/combigrid/grid/CombinationGrid.hpp>
#include <sgpp/combigrid/grid/IndexVectorRange.hpp>

#include <algorithm>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sgpp {
namespace combigrid {

namespace {
class HashLevelVector {
 public:
  size_t operator()(const LevelVector& level) const {
    size_t seed = level.size();

    for (const level_t& l : level) {
      seed ^= (l + 0x9e3779b9 + (seed << 6) + (seed >> 2));
    }

    return seed;
  }
};
}  // namespace

CombinationGrid::CombinationGrid() : fullGrids(), coefficients() {
}

CombinationGrid::CombinationGrid(const std::vector<FullGrid>& fullGrids,
    const base::DataVector& coefficients) : fullGrids(fullGrids), coefficients(coefficients) {
}

CombinationGrid CombinationGrid::fromRegular(size_t dim, level_t n,
    const HeterogeneousBasis& basis, bool hasBoundary) {
  std::vector<size_t> binomialCoefficients((dim+1)/2);
  binomialCoefficients[0] = 1.0;

  // binomial(dim-1, d)
  // = ((dim-1) * (dim-2) * ... * (dim-d)) / (1 * 2 * ... * d)
  // = binomial(dim-1, d-1) * (dim-d) / d
  for (size_t q = 1; q < (dim+1)/2; q++) {
    binomialCoefficients[q] = binomialCoefficients[q-1] * (dim-q) / q;
  }

  std::vector<FullGrid> fullGrids;
  base::DataVector coefficients;
  const level_t maxLevelSum = (hasBoundary ? n : static_cast<level_t>(n+dim-1));

  for (size_t q = 0; q < dim; q++) {
    const std::vector<LevelVector> levels = (hasBoundary ?
        enumerateLevelsWithSumWithBoundary(dim, maxLevelSum-static_cast<level_t>(q)) :
        enumerateLevelsWithSumWithoutBoundary(dim, maxLevelSum-static_cast<level_t>(q)));
    const double coefficient = ((q % 2 == 0) ? 1.0 : -1.0) *
        static_cast<double>(binomialCoefficients[((q < (dim+1)/2) ? q : (dim-q-1))]);

    for (const LevelVector& level : levels) {
      fullGrids.emplace_back(level, basis, hasBoundary);
      coefficients.push_back(coefficient);
    }
  }

  return CombinationGrid(fullGrids, coefficients);
}

CombinationGrid CombinationGrid::fromSubspaces(
    const std::vector<LevelVector>& subspaceLevels, const HeterogeneousBasis& basis,
    bool hasBoundary) {
  const size_t dim = basis.getDimension();
  std::unordered_map<LevelVector, double, HashLevelVector> coefficientMap;
  LevelVector offsetMinIndex(dim, 0);
  LevelVector offsetMaxIndex(dim, 1);
  LevelVector level(dim);
  std::vector<FullGrid> fullGrids;
  base::DataVector coefficients;

  for (const LevelVector& subspaceLevel : subspaceLevels) {
    IndexVectorRange indexVectorRange(offsetMinIndex, offsetMaxIndex);
    double coefficient = 0.0;

    for (const IndexVector& offsetIndex : indexVectorRange) {
      double sign = 1.0;

      for (size_t d = 0; d < dim; d++) {
        level[d] = subspaceLevel[d] + offsetIndex[d];
        sign *= (1.0 - 2.0 * static_cast<double>(offsetIndex[d]));
      }

      if (std::find(subspaceLevels.begin(), subspaceLevels.end(), level) != subspaceLevels.end()) {
        coefficient += sign;
      }
    }

    if (coefficient != 0.0) {
      fullGrids.emplace_back(subspaceLevel, basis, hasBoundary);
      coefficients.push_back(coefficient);
    }
  }

  return CombinationGrid(fullGrids, coefficients);
}

void CombinationGrid::combinePoints(base::GridStorage& gridStorage) const {
  const size_t dim = getDimension();
  base::GridPoint point(dim);
  gridStorage.clear();

  for (const FullGrid& fullGrid : fullGrids) {
    const LevelVector& level = fullGrid.getLevel();

    for (const IndexVector& index : IndexVectorRange(fullGrid)) {
      for (size_t d = 0; d < dim; d++) {
        level_t l = level[d];
        index_t i = index[d];
        HeterogeneousBasis::hierarchizeLevelIndex(l, i);
        point.set(d, l, i);
      }

      if (!gridStorage.isContaining(point)) {
        gridStorage.insert(point);
      }
    }
  }
}

double CombinationGrid::combineValues(const base::DataVector& values) const {
  return values.dotProduct(coefficients);
}

void CombinationGrid::combineValues(const base::DataMatrix& values,
    base::DataVector& result) const {
  result.resize(values.getNrows());
  values.mult(coefficients, result);
}

void CombinationGrid::combineSparseGridValues(const base::GridStorage& gridStorage,
    const std::vector<base::DataVector>& values, base::DataVector& result) const {
  const size_t N = gridStorage.getSize();
  const size_t n = fullGrids.size();
  const size_t dim = getDimension();
  IndexVector index(dim);
  IndexVectorRange range;
  result.resize(N);
  result.setAll(0.0);

  for (size_t k = 0; k < N; k++) {
    for (size_t i = 0; i < n; i++) {
      if (findGridPointInFullGrid(fullGrids[i], gridStorage[k], index)) {
        range.setGrid(fullGrids[i]);
        result[k] += coefficients[i] * values[i][range.find(index)];
      }
    }
  }
}

void CombinationGrid::combineSparseGridValues(const base::GridStorage& gridStorage,
    const std::vector<base::DataMatrix>& values, base::DataMatrix& result) const {
  const size_t N = gridStorage.getSize();

  if (values.empty()) {
    result.resize(N, 0);
    return;
  }

  const size_t n = fullGrids.size();
  const size_t m = values[0].getNrows();
  const size_t dim = getDimension();
  IndexVector index(dim);
  IndexVectorRange range;
  result.resize(N, m);
  result.setAll(0.0);

  for (size_t k = 0; k < N; k++) {
    for (size_t i = 0; i < n; i++) {
      if (findGridPointInFullGrid(fullGrids[i], gridStorage[k], index)) {
        range.setGrid(fullGrids[i]);

        for (size_t j = 0; j < m; j++) {
          result(k, j) += coefficients[i] * values[i](range.find(index), j);
        }
      }
    }
  }
}

void CombinationGrid::distributeValuesToFullGrid(const base::GridStorage& gridStorage,
    const base::DataVector& values, const FullGrid& fullGrid, base::DataVector& result) const {
  const size_t dim = getDimension();
  const size_t N = gridStorage.getSize();
  IndexVector index(dim);
  IndexVectorRange range(fullGrid);
  result.resize(fullGrid.getNumberOfIndexVectors());
  result.setAll(0.0);

  for (size_t k = 0; k < N; k++) {
    if (findGridPointInFullGrid(fullGrid, gridStorage[k], index)) {
      result[range.find(index)] = values[k];
    }
  }
}

void CombinationGrid::distributeValuesToFullGrids(const base::GridStorage& gridStorage,
    const base::DataVector& values, std::vector<base::DataVector>& result) const {
  const size_t dim = getDimension();
  const size_t N = gridStorage.getSize();
  const size_t n = fullGrids.size();
  IndexVector index(dim);
  IndexVectorRange range;
  result.resize(n);

  for (size_t i = 0; i < n; i++) {
    result[i].resize(fullGrids[i].getNumberOfIndexVectors());
    result[i].setAll(0.0);
  }

  for (size_t k = 0; k < N; k++) {
    for (size_t i = 0; i < n; i++) {
      if (findGridPointInFullGrid(fullGrids[i], gridStorage[k], index)) {
        range.setGrid(fullGrids[i]);
        result[i][range.find(index)] = values[k];
      }
    }
  }
}

size_t CombinationGrid::getDimension() const {
  return (fullGrids.empty() ? 0 : fullGrids[0].getDimension());
}

const std::vector<FullGrid>& CombinationGrid::getFullGrids() const {
  return fullGrids;
}

const base::DataVector& CombinationGrid::getCoefficients() const {
  return coefficients;
}

void CombinationGrid::setFullGridsAndCoefficients(const std::vector<FullGrid>& fullGrids,
    const base::DataVector& coefficients) {
  this->fullGrids = fullGrids;
  this->coefficients = coefficients;
}

std::vector<LevelVector> CombinationGrid::enumerateLevelsWithSumWithBoundary(
    size_t dim, level_t n) {
  if (dim == 0) {
    return std::vector<LevelVector>();
  } else if (dim == 1) {
    return std::vector<LevelVector>{LevelVector{n}};
  } else if (n == 0) {
    return std::vector<LevelVector>{LevelVector(dim, 0)};
  } else {
    std::vector<LevelVector> result;

    for (level_t l = 0; l <= n; l++) {
      for (LevelVector& level : enumerateLevelsWithSumWithBoundary(dim-1, n-l)) {
        level.push_back(l);
        result.push_back(level);
      }
    }

    return result;
  }
}

std::vector<LevelVector> CombinationGrid::enumerateLevelsWithSumWithoutBoundary(
    size_t dim, level_t n) {
  if ((dim == 0) || (n < dim)) {
    return std::vector<LevelVector>();
  } else if (dim == 1) {
    return std::vector<LevelVector>{LevelVector{n}};
  } else {
    std::vector<LevelVector> result;

    for (level_t l = 1; l <= n-dim+1; l++) {
      for (LevelVector& level : enumerateLevelsWithSumWithoutBoundary(dim-1, n-l)) {
        level.push_back(l);
        result.push_back(level);
      }
    }

    return result;
  }
}

bool CombinationGrid::findGridPointInFullGrid(const FullGrid& fullGrid,
    const base::GridPoint& gridPoint, IndexVector& index) {
  const LevelVector& levelFullGrid = fullGrid.getLevel();
  bool isContained = true;

  for (size_t d = 0; d < gridPoint.getDimension(); d++) {
    const level_t curLevel = gridPoint.getLevel(d);

    if ((curLevel <= levelFullGrid[d]) && (fullGrid.hasBoundary() || (curLevel >= 1))) {
      index[d] = gridPoint.getIndex(d) << (levelFullGrid[d] - curLevel);
    } else {
      isContained = false;
      break;
    }
  }

  return isContained;
}

}  // namespace combigrid
}  // namespace sgpp
