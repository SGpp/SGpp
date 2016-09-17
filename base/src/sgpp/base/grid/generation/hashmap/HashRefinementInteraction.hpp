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

class HashRefinementInteraction : public HashRefinement {
 public:
  explicit HashRefinementInteraction(std::unordered_set<std::vector<bool>> interactions);
  void createGridpoint(GridStorage& storage, index_type& index) override;

 private:
  std::unordered_set<std::vector<bool>> interactions;
};

}  // namespace base
}  // namespace sgpp

#endif  // HASHREFINEMENTINTERACTION_HPP
