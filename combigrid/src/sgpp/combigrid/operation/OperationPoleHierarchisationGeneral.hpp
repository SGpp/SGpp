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
#include <sgpp/combigrid/operation/OperationPole.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class OperationPoleHierarchisationGeneral : public OperationPole {
 public:
  explicit OperationPoleHierarchisationGeneral(base::Basis<level_t, index_t>& basis,
      bool isBasisHierarchical);

  ~OperationPoleHierarchisationGeneral() override;

  static void fromHeterogenerousBasis(const HeterogeneousBasis& basis,
      std::vector<std::unique_ptr<OperationPole>>& operation);

  static void fromHeterogenerousBasis(const HeterogeneousBasis& basis,
      std::vector<OperationPole*>& operation);

  void apply(base::DataVector& values, size_t start, size_t step, size_t count,
      level_t level, bool hasBoundary = true) override;

 protected:
  class HierarchisationGeneralSLE : public base::SLE {
   public:
    HierarchisationGeneralSLE(base::Basis<level_t, index_t>& basis, bool isBasisHierarchical,
        size_t dim, level_t level, bool hasBoundary = true);

    ~HierarchisationGeneralSLE() override;

    double getMatrixEntry(size_t i, size_t j) override;

    bool isMatrixEntryNonZero(size_t i, size_t j) override;

    bool isBasisHierarchical() const;

    void setIsBasisHierarchical(bool isBasisHierarchical);

    size_t getDimension() const override;

    void setDimension(size_t dim);

    level_t getLevel() const;

    void setLevel(level_t level);

    bool hasBoundary() const;

    void setHasBoundary(bool hasBoundary);

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
