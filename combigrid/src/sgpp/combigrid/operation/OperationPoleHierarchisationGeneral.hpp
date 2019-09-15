// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/tools/sle/solver/Auto.hpp>
#include <sgpp/base/tools/sle/system/SLE.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/basis/HeterogeneousBasis.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>
#include <sgpp/combigrid/operation/OperationPole.hpp>

#include <cmath>
#include <vector>

namespace sgpp {
namespace combigrid {

class OperationPoleHierarchisationGeneral : public OperationPole {
 public:
  explicit OperationPoleHierarchisationGeneral(base::Basis<level_t, index_t>& basis,
      bool isBasisHierarchical) : sle(basis, isBasisHierarchical, 0, 0), sleSolver() {
  }

  ~OperationPoleHierarchisationGeneral() override {
  }

  static void fromHeterogenerousBasis(const HeterogeneousBasis& basis,
      std::vector<std::unique_ptr<OperationPole>>& operation) {
    for (base::Basis<level_t, index_t>* const & basis1d : basis.getBases1d()) {
      operation.emplace_back(new OperationPoleHierarchisationGeneral(
          *basis1d, basis.isHierarchical()));
    }
  }

  static void fromHeterogenerousBasis(const HeterogeneousBasis& basis,
      std::vector<OperationPole*>& operation) {
    for (base::Basis<level_t, index_t>* const & basis1d : basis.getBases1d()) {
      operation.push_back(new OperationPoleHierarchisationGeneral(
          *basis1d, basis.isHierarchical()));
    }
  }

  void apply(base::DataVector& values, size_t start, size_t step, size_t count,
      level_t level, bool hasBoundary = true) override {
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

 protected:
  class HierarchisationGeneralSLE : public base::SLE {
   public:
    HierarchisationGeneralSLE(base::Basis<level_t, index_t>& basis, bool isBasisHierarchical,
        size_t dim, level_t level, bool hasBoundary = true) :
        basis(basis), isBasisHierarchical_(isBasisHierarchical),
        dim(dim), level(level), hasBoundary_(hasBoundary) {
    }

    ~HierarchisationGeneralSLE() override {
    }

    double getMatrixEntry(size_t i, size_t j) override {
      level_t levelBasis = level;
      index_t indexBasis = static_cast<index_t>(j);

      if (isBasisHierarchical_) {
        HeterogeneousBasis::hierarchizeLevelIndex(levelBasis, indexBasis);
      }

      const double point = static_cast<double>(i + (hasBoundary_ ? 0 : 1)) /
          (static_cast<index_t>(1) << level);
      return basis.eval(levelBasis, indexBasis, point);
    }

    bool isMatrixEntryNonZero(size_t i, size_t j) override {
      return (std::abs(getMatrixEntry(i, j)) > 1e-12);
    }

    bool isBasisHierarchical() const {
      return isBasisHierarchical_;
    }

    void setIsBasisHierarchical(bool isBasisHierarchical) {
      isBasisHierarchical_ = isBasisHierarchical;
    }

    size_t getDimension() const override {
      return dim;
    }

    void setDimension(size_t dim) {
      this->dim = dim;
    }

    level_t getLevel() const {
      return level;
    }

    void setLevel(level_t level) {
      this->level = level;
    }

    bool hasBoundary() const {
      return hasBoundary_;
    }

    void setHasBoundary(bool hasBoundary) {
      hasBoundary_ = hasBoundary;
    }

   protected:
    base::Basis<level_t, index_t>& basis;
    bool isBasisHierarchical_;
    size_t dim;
    level_t level;
    bool hasBoundary_;
  };

  HierarchisationGeneralSLE sle;
  base::sle_solver::Auto sleSolver;
};

}  // namespace combigrid
}  // namespace sgpp
