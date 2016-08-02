// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HASHREFINEMENTINCONSISTENT_HPP
#define HASHREFINEMENTINCONSISTENT_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Free refinement class for sparse grids
 */
class HashRefinementInconsistent: public HashRefinement {
 protected:
  /**
   * This method creates a new point on the grid. It checks if some parents or
   * children are needed in other dimensions.
   *
   * @param storage hashmap that stores the gridpoints
   * @param point The point that should be inserted
   */
  void createGridpoint(GridStorage& storage, GridPoint& point) override;
};

}  // namespace base
}  // namespace sgpp

#endif /* HASHREFINEMENTINCONSISTENT_HPP */
