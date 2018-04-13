// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ANOVAREFINEMENT_HPP_
#define ANOVAREFINEMENT_HPP_

#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/*
 * ANOVAHashRefinement inserts creates children only in the dimensions where the
 * level is greater than 1.
 */
class ANOVAHashRefinement : public HashRefinement {
 public:
  /**
     * This method refines a grid point by generating the children in every dimension
     * of the grid where the level is greater than 1 and all their missing ancestors by calling
     * create_gridpoint().
     *
     * @param storage hashmap that stores the gridpoints
     * @param refine_index The index in the hashmap of the point that should be refined
     */
  virtual void refineGridpoint(GridStorage& storage, size_t refine_index);
};
}  // namespace base
}  // namespace sgpp

#endif /* ANOVAREFINEMENT_HPP_ */
