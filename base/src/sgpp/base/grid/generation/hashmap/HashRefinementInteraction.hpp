// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HASHREFINEMENTINTERACTION_HPP
#define HASHREFINEMENTINTERACTION_HPP

#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <unordered_set>
#include <vector>

namespace sgpp {
namespace base {
/**
 * @brief The HashRefinementInteraction class
 * @details Refines the grid, but only adds interactions that are contained in the set interactions, i.e.
 * only desired interactions.
 * Each desired interaction is encoded as a vector, of which each entry is true if the dimension
 * should be used and false otherwise.
 * For example, a \f$ (x_1 \times x_2)\f$ interaction for a 3-dimensional dataset is
 * modelled as the boolean vector (true, true, false).
 */
class HashRefinementInteraction : public HashRefinement {
 public:
    /**
   * @brief HashRefinementInteraction
   * @param interactions contains all desired interactions.
   */
  explicit HashRefinementInteraction(std::unordered_set<std::vector<bool>> interactions);
  void createGridpoint(GridStorage& storage, index_type& index) override;

 private:
  std::unordered_set<std::vector<bool>> interactions;
};

}  // namespace base
}  // namespace sgpp

#endif  // HASHREFINEMENTINTERACTION_HPP
