// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/hashmap/HashRefinementMultipleClass.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/tools/MultipleClassPoint.hpp>

#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include <cmath>

namespace sgpp {
namespace base {

HashRefinementMultipleClass::HashRefinementMultipleClass(Grid& grid,
        std::vector<sgpp::base::MultipleClassPoint>* pts,
        std::vector<Grid*>& classGrids, double &borderSum, double &borderCnt,
        double topPercent) : HashRefinement(), points(pts),
        multigrid(grid), grids(classGrids), borderSum(borderSum),
        borderCnt(borderCnt), topPercent(topPercent) {
}

void HashRefinementMultipleClass::refineGridpoint(GridStorage& storage,
                                     size_t refine_index) {
    // find index in combined grid
    GridPoint point(storage[refine_index]);
    size_t multiSeq = multigrid.getStorage().getSequenceNumber(point);

    if ( multiSeq >= points->size() ) {
        // new point, not in combined grid
        return;
    }
    sgpp::base::MultipleClassPoint mcp = points->at(multiSeq);
    GridStorage& tStorage = grids.at(mcp.getDominateClass())->getStorage();

    // add top classes
    std::vector<std::tuple<double, size_t, bool>> top = mcp.getTopClasses(topPercent);
    for (size_t n = 0 ; n < top.size() ; n++) {
        // add points to all top classes
        GridStorage& nStorage = grids.at(std::get<1>(top.at(n)))->getStorage();
        addGridpoint(nStorage, point);
    }

    // add neighbors
    std::vector<std::tuple<size_t, size_t, bool>> neighbors = mcp.getNeighbors();
    for (size_t n = 0 ; n < neighbors.size() ; n++) {
        size_t nextPoint = std::get<0>(neighbors.at(n));
        sgpp::base::MultipleClassPoint neighP = points->at(nextPoint);
        size_t dim = std::get<1>(neighbors.at(n));
        bool isLeft = std::get<2>(neighbors.at(n));

        GridStorage& nStorage = grids.at(neighP.getDominateClass())->getStorage();
        GridPoint nPoint = multigrid.getStorage().getPoint(nextPoint);

        // adds original point to needed grid
        addGridpoint(tStorage, point);
        addGridpoint(nStorage, point);

        // adds neighbor to needed grid
        addGridpoint(tStorage, nPoint);
        addGridpoint(nStorage, nPoint);

        // create children in given dimension/direction
        // from refineGridpoint1D
        index_t source_index;
        level_t source_level;
        if ( point.getLevel(dim) < nPoint.getLevel(dim) ) {
            nPoint.get(dim, source_level, source_index);
            isLeft = !isLeft;
        } else {
            point.get(dim, source_level, source_index);
        }
        if ( isLeft ) {
            // generate left child, if necessary
            point.set(dim, source_level + 1, 2 * source_index - 1);
        } else {
            // generate right child, if necessary
            point.set(dim, source_level + 1, 2 * source_index + 1);
        }
        // insert point in both dominate classes
        addGridpoint(tStorage, point);
        addGridpoint(nStorage, point);

        point.set(dim, source_level, source_index);
    }

    // add points to border
    std::vector<std::tuple<size_t, size_t, bool>> borders = mcp.getBorders();
    if ( borderSum/borderCnt <
            mcp.getBorderScore() * mcp.getDensity(mcp.getDominateClass())) {
        for (size_t b = 0 ; b < borders.size() ; b++) {
            size_t dim = std::get<0>(borders.at(b));
            bool isLeft = std::get<2>(borders.at(b));

            index_t source_index;
            level_t source_level;
            point.get(dim, source_level, source_index);

            if ( isLeft ) {
                // generate left child, if necessary
                point.set(dim, source_level + 1, 2 * source_index - 1);
            } else {
                // generate right child, if necessary
                point.set(dim, source_level + 1, 2 * source_index + 1);
            }
            addGridpoint(grids.at(mcp.getDominateClass())->getStorage(), point);

            point.set(dim, source_level, source_index);
        }
    }
}

void HashRefinementMultipleClass::addGridpoint(GridStorage& storage, GridPoint& point) {
    if (!storage.isContaining(point)) {
        createGridpoint(storage, point);
    }
}

void HashRefinementMultipleClass::collectRefinablePoints(GridStorage& storage,
    RefinementFunctor& functor,
    AbstractRefinement::refinement_container_type& collection) {
  size_t refinements_num = functor.getRefinementsNum();

  // max value equals min value
  GridPoint point;
  GridStorage::grid_map_iterator end_iter = storage.end();

  // start iterating over whole grid
  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter;
       iter++) {
    point = *(iter->first);
    GridStorage::grid_map_iterator child_iter;

    // all points can be refined
    AbstractRefinement::refinement_list_type current_value_list =
        getIndicator(storage, iter, functor);
    addElementToCollection(iter, current_value_list, refinements_num,
        collection);
  }
}

} /* namespace base */
} /* namespace sgpp */
