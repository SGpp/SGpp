// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/combigrid/grid/CombinationGrid.hpp>
#include <sgpp/combigrid/tools/IndexVectorRange.hpp>
#include <sgpp/combigrid/tools/LevelVectorTools.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

namespace sgpp {
namespace combigrid {

CombinationGrid::CombinationGrid() : fullGrids(), coefficients() {}

CombinationGrid::CombinationGrid(const std::vector<FullGrid>& fullGrids,
                                 const base::DataVector& coefficients)
    : fullGrids(fullGrids), coefficients(coefficients) {}

CombinationGrid::CombinationGrid(const FullGrid& fullGrid)
    : fullGrids(std::vector<FullGrid>(1, fullGrid)), coefficients(base::DataVector(1, 1)) {}

CombinationGrid CombinationGrid::fromRegularSparse(size_t dim, level_t n,
                                                   const HeterogeneousBasis& basis,
                                                   bool hasBoundary) {
  return fromRegularSparseTruncated(dim, LevelVector(dim, 0), n, basis, hasBoundary);
}

CombinationGrid CombinationGrid::fromRegularSparseTruncated(size_t dim, LevelVector truncationLevel,
                                                            level_t levelSumDistance,
                                                            const HeterogeneousBasis& basis,
                                                            bool hasBoundary) {
  assert(truncationLevel.size() == dim);

  std::vector<size_t> binomialCoefficients((dim + 1) / 2);
  binomialCoefficients[0] = 1.0;

  // binomial(dim-1, d)
  // = ((dim-1) * (dim-2) * ... * (dim-d)) / (1 * 2 * ... * d)
  // = binomial(dim-1, d-1) * (dim-d) / d
  for (size_t q = 1; q < (dim + 1) / 2; q++) {
    binomialCoefficients[q] = binomialCoefficients[q - 1] * (dim - q) / q;
  }

  std::vector<FullGrid> fullGrids;
  base::DataVector coefficients;
  const level_t maxLevelDifference =
      (hasBoundary ? levelSumDistance : static_cast<level_t>(levelSumDistance + dim - 1));

  for (size_t q = 0; q < dim; q++) {
    if (q <= maxLevelDifference) {
      std::vector<LevelVector> levels =
          (hasBoundary ? LevelVectorTools::generateDiagonalWithBoundary(
                             dim, maxLevelDifference - static_cast<level_t>(q))
                       : LevelVectorTools::generateDiagonalWithoutBoundary(
                             dim, maxLevelDifference - static_cast<level_t>(q)));
      const double coefficient =
          ((q % 2 == 0) ? 1.0 : -1.0) *
          static_cast<double>(binomialCoefficients[((q < (dim + 1) / 2) ? q : (dim - q - 1))]);

      for (LevelVector& level : levels) {
        // difference between truncated and regular combination technique: offset by the minimum
        // level / truncation level
        for (size_t d = 0; d < dim; ++d) {
          level[d] += truncationLevel[d];
        }
        fullGrids.emplace_back(level, basis, hasBoundary);
        coefficients.push_back(coefficient);
      }
    }
  }

  return CombinationGrid(fullGrids, coefficients);
}

CombinationGrid CombinationGrid::fromSubspaces(const std::vector<LevelVector>& subspaceLevels,
                                               const HeterogeneousBasis& basis, bool hasBoundary) {
  std::vector<FullGrid> fullGrids;
  auto coefficients = getStandardCoefficientsFromLevelSet(subspaceLevels);
  decltype(coefficients) coefficientsNonZero;

  for (size_t i = 0; i < subspaceLevels.size(); ++i) {
    if (coefficients[i] != 0.0) {
      fullGrids.emplace_back(subspaceLevels[i], basis, hasBoundary);
      coefficientsNonZero.push_back(coefficients[i]);
    }
  }
  return CombinationGrid(fullGrids, coefficientsNonZero);
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
                                              const std::vector<base::DataVector>& values,
                                              base::DataVector& result) const {
  const size_t N = gridStorage.getSize();
  const size_t n = fullGrids.size();
  const size_t dim = getDimension();
  IndexVector index(dim);
  IndexVectorRange range;
  result.resize(N);
  result.setAll(0.0);

  for (size_t k = 0; k < N; k++) {
    for (size_t i = 0; i < n; i++) {
      if (fullGrids[i].findGridPointInFullGrid(gridStorage[k], index)) {
        range.setGrid(fullGrids[i]);
        result[k] += coefficients[i] * values[i][range.find(index)];
      }
    }
  }
}

void CombinationGrid::combineSparseGridValues(const base::GridStorage& gridStorage,
                                              const std::vector<base::DataMatrix>& values,
                                              base::DataMatrix& result) const {
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
      if (fullGrids[i].findGridPointInFullGrid(gridStorage[k], index)) {
        range.setGrid(fullGrids[i]);

        for (size_t j = 0; j < m; j++) {
          result(k, j) += coefficients[i] * values[i](range.find(index), j);
        }
      }
    }
  }
}

void CombinationGrid::distributeValuesToFullGrid(const base::GridStorage& gridStorage,
                                                 const base::DataVector& values,
                                                 const FullGrid& fullGrid,
                                                 base::DataVector& result) const {
  const size_t dim = getDimension();
  const size_t N = gridStorage.getSize();
  IndexVector index(dim);
  IndexVectorRange range(fullGrid);
  result.resize(fullGrid.getNumberOfIndexVectors());
  result.setAll(0.0);

  for (size_t k = 0; k < N; k++) {
    if (fullGrid.findGridPointInFullGrid(gridStorage[k], index)) {
      result[range.find(index)] = values[k];
    }
  }
}

void CombinationGrid::distributeValuesToFullGrids(const base::GridStorage& gridStorage,
                                                  const base::DataVector& values,
                                                  std::vector<base::DataVector>& result) const {
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
      if (fullGrids[i].findGridPointInFullGrid(gridStorage[k], index)) {
        range.setGrid(fullGrids[i]);
        result[i][range.find(index)] = values[k];
      }
    }
  }
}

size_t CombinationGrid::getDimension() const {
  return (fullGrids.empty() ? 0 : fullGrids[0].getDimension());
}

const std::vector<FullGrid>& CombinationGrid::getFullGrids() const { return fullGrids; }

const base::DataVector& CombinationGrid::getCoefficients() const { return coefficients; }

void CombinationGrid::setFullGridsAndCoefficients(const std::vector<FullGrid>& fullGrids,
                                                  const base::DataVector& coefficients) {
  this->fullGrids = fullGrids;
  this->coefficients = coefficients;
}

base::DataVector getStandardCoefficientsFromLevelSet(const std::vector<LevelVector>& levelSet) {
  const size_t dim = levelSet[0].size();
  LevelVector offsetMinIndex(dim, 0);
  LevelVector offsetMaxIndex(dim, 1);
  LevelVector tmp_level(dim);
  base::DataVector coefficients;
  coefficients.reserve(levelSet.size());

  for (const LevelVector& level : levelSet) {
    IndexVectorRange indexVectorRange(offsetMinIndex, offsetMaxIndex);
    double coefficient = 0.0;

    for (const IndexVector& offsetIndex : indexVectorRange) {
      double sign = 1.0;

      for (size_t d = 0; d < dim; d++) {
        tmp_level[d] = level[d] + offsetIndex[d];
        sign *= (1.0 - 2.0 * static_cast<double>(offsetIndex[d]));
      }

      if (std::find(levelSet.begin(), levelSet.end(), tmp_level) != levelSet.end()) {
        coefficient += sign;
      }
    }
    coefficients.push_back(coefficient);
  }
  return coefficients;
}

}  // namespace combigrid
}  // namespace sgpp
