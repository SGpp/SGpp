// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HASHREFINEMENTMULTIPLECLASS_HPP
#define HASHREFINEMENTMULTIPLECLASS_HPP

#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/tools/MultipleClassPoint.hpp>

#include <iostream>
#include <tuple>
#include <cmath>
#include <vector>
#include <algorithm>


namespace sgpp {
namespace base {

class HashRefinementMultipleClass : public HashRefinement {
 public:
    HashRefinementMultipleClass(Grid& grid,
        std::vector<sgpp::base::MultipleClassPoint>* pts,
        std::vector<Grid*>& classGrids,
        double &borderSum, double &borderCnt, double topPercent);
    virtual ~HashRefinementMultipleClass() {}

 protected:
  void refineGridpointsCollection(
        GridStorage& storage,
        RefinementFunctor& functor,
        AbstractRefinement::refinement_container_type& collection) override;
  void refineGridpoint(GridStorage& storage, size_t refine_index) override;
  void collectRefinablePoints(GridStorage& storage,
        RefinementFunctor& functor,
        AbstractRefinement::refinement_container_type& collection) override;

 private:
    std::vector<sgpp::base::MultipleClassPoint>* points;
    Grid& multigrid;
    std::vector<Grid*>& grids;
    double &borderSum;
    double &borderCnt;
    double topPercent;

    void addGridpoint(GridStorage& storage, GridPoint& point);
};

} /* namespace base */
} /* namespace sgpp */

#endif /* HASHREFINEMENTMULTIPLECLASS_HPP */
