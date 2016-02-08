// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SQUAREROOTGRIDGENERATOR_HPP_
#define SQUAREROOTGRIDGENERATOR_HPP_
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * This class provides the interface for the grid generation
 * for grids with boundaries, pentagon cut through sub space scheme
 */
class SquareRootGridGenerator : public GridGenerator {
 public:
  /**
   * Constructor
   *
   * @param storage template type that holds the grid points
   */
  SquareRootGridGenerator(GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~SquareRootGridGenerator();

  virtual void regular(size_t level) override;
  virtual void cliques(size_t level, size_t clique_size) override;
  virtual void full(size_t level) override {};
  virtual void refine(RefinementFunctor* func) override {};
  virtual size_t getNumberOfRefinablePoints() override {
    return 0;
  };

  virtual void coarsen(CoarseningFunctor* func, DataVector* alpha) override {};
  virtual void coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha,
                                 size_t numFirstOnly) override {};
  virtual size_t getNumberOfRemovablePoints() override {
    return 0;
  };

  virtual void refineMaxLevel(RefinementFunctor* func,
                              size_t maxLevel) override {};
  virtual size_t getNumberOfRefinablePointsToMaxLevel(size_t maxLevel) override {
    return 0;
  };

 protected:
  /// Pointer to the grid's storage object
  GridStorage* storage;
};

}  // namespace base
}  // namespace SGPP

#endif /* SQUAREROOTGRIDGENERATOR_HPP_ */
