// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/combigrid/operation/OperationPoleHierarchisationGeneral.hpp>

#include <cmath>
#include <vector>

namespace sgpp {
namespace combigrid {

OperationPoleHierarchisationGeneral::OperationPoleHierarchisationGeneral(
    base::Basis<level_t, index_t>& basis, bool isBasisHierarchical) :
    sle(basis, 0, 0, isBasisHierarchical), sleSolver() {
}

OperationPoleHierarchisationGeneral::~OperationPoleHierarchisationGeneral() {
}

void OperationPoleHierarchisationGeneral::fromHeterogenerousBasis(const HeterogeneousBasis& basis,
    std::vector<std::unique_ptr<OperationPole>>& operation) {
  for (base::Basis<level_t, index_t>* const & basis1d : basis.getBases1d()) {
    operation.emplace_back(new OperationPoleHierarchisationGeneral(
        *basis1d, basis.isHierarchical()));
  }
}

void OperationPoleHierarchisationGeneral::fromHeterogenerousBasis(const HeterogeneousBasis& basis,
    std::vector<OperationPole*>& operation) {
  for (base::Basis<level_t, index_t>* const & basis1d : basis.getBases1d()) {
    operation.push_back(new OperationPoleHierarchisationGeneral(
        *basis1d, basis.isHierarchical()));
  }
}

void OperationPoleHierarchisationGeneral::apply(base::DataVector& values, size_t start, size_t step,
    size_t count, level_t level, bool hasBoundary) {
  base::DataVector rhs(count);
  base::DataVector solution(count);

  for (size_t i = 0; i < count; i++) {
    rhs[i] = values[start + i * step];
  }

  sle.setDimension(count);
  sle.setLevel(level);
  sle.setHasBoundary(hasBoundary);
  sleSolver.solve(sle, rhs, solution);

  for (size_t i = 0; i < count; i++) {
    values[start + i * step] = solution[i];
  }
}

OperationPoleHierarchisationGeneral::HierarchisationGeneralSLE::HierarchisationGeneralSLE(
    base::Basis<level_t, index_t>& basis, size_t dim, level_t level,
    bool isBasisHierarchical, bool hasBoundary) :
    basis(basis), isBasisHierarchical_(isBasisHierarchical),
    dim(dim), level(level), hasBoundary_(hasBoundary) {
}

OperationPoleHierarchisationGeneral::HierarchisationGeneralSLE::~HierarchisationGeneralSLE() {
}

double OperationPoleHierarchisationGeneral::HierarchisationGeneralSLE::getMatrixEntry(
    size_t i, size_t j) {
  level_t levelBasis = level;
  index_t indexBasis = static_cast<index_t>(j);

  if (isBasisHierarchical_) {
    HeterogeneousBasis::hierarchizeLevelIndex(levelBasis, indexBasis);
  }

  const double point = static_cast<double>(i + (hasBoundary_ ? 0 : 1)) /
      (static_cast<index_t>(1) << level);
  return basis.eval(levelBasis, indexBasis, point);
}

bool OperationPoleHierarchisationGeneral::HierarchisationGeneralSLE::isMatrixEntryNonZero(
    size_t i, size_t j) {
  return (std::abs(getMatrixEntry(i, j)) > 1e-12);
}

bool OperationPoleHierarchisationGeneral::HierarchisationGeneralSLE::isBasisHierarchical() const {
  return isBasisHierarchical_;
}

void OperationPoleHierarchisationGeneral::HierarchisationGeneralSLE::setIsBasisHierarchical(
    bool isBasisHierarchical) {
  isBasisHierarchical_ = isBasisHierarchical;
}

size_t OperationPoleHierarchisationGeneral::HierarchisationGeneralSLE::getDimension() const {
  return dim;
}

void OperationPoleHierarchisationGeneral::HierarchisationGeneralSLE::setDimension(size_t dim) {
  this->dim = dim;
}

level_t OperationPoleHierarchisationGeneral::HierarchisationGeneralSLE::getLevel() const {
  return level;
}

void OperationPoleHierarchisationGeneral::HierarchisationGeneralSLE::setLevel(level_t level) {
  this->level = level;
}

bool OperationPoleHierarchisationGeneral::HierarchisationGeneralSLE::hasBoundary() const {
  return hasBoundary_;
}

void OperationPoleHierarchisationGeneral::HierarchisationGeneralSLE::setHasBoundary(
    bool hasBoundary) {
  hasBoundary_ = hasBoundary;
}

}  // namespace combigrid
}  // namespace sgpp
