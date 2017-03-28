/*
 * HashRefinementMultipleClass.cpp
 *
 *  Created on: Mar 9, 2017
 *      Author: katrin
 */

#include "HashRefinementMultipleClass.hpp"
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <iostream>
#include <tuple>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sgpp/base/tools/MultipleClassPoint.hpp>

namespace sgpp {
namespace base {

    // TODO (degel_kn): general
HashRefinementMultipleClass::HashRefinementMultipleClass(Grid& grid,
        std::vector<sgpp::base::MultipleClassPoint>& pts,
        std::vector<Grid*>& classGrids) : HashRefinement(),
        points(pts), multigrid(grid), grids(classGrids) {
    
    /*grids = classGrids;
    multigrid = grid;
    points = pts;*/
}

void HashRefinementMultipleClass::refineGridpoint(GridStorage& storage,
                                     size_t refine_index) {    
    // find index in combined grid
    std::cout << " refineGridpoint: " << refine_index << std::endl;
    GridPoint point(storage[refine_index]);
    size_t multiSeq = multigrid.getStorage().getSequenceNumber(point);

    if ( multiSeq >= points.size() ) {
        // new point, not in combined grid
        return;
    }
    std::string pInfo1 = "(";
    for (size_t d = 0; d < multigrid.getStorage().getDimension(); d++) {
        pInfo1 += std::to_string(point.getStandardCoordinate(d)) + ",";
    }
    GridPoint xPoint = multigrid.getStorage().getPoint(multiSeq);      
    std::string pInfo2 = "(";
    for (size_t d = 0; d < multigrid.getStorage().getDimension(); d++) {
        pInfo2 += std::to_string(xPoint.getStandardCoordinate(d)) + ",";
    }
    std::cout << "refine point: " << refine_index << " -- " << pInfo1
            << ") << (" << multiSeq << "/" << points.size() << ")" << " -- " << pInfo2 << std::endl;

    // GridPoint point(storage[refine_index]);
    // Sets leaf property of index, which is refined to false
    // storage[refine_index].setLeaf(false);

    sgpp::base::MultipleClassPoint mcp = points.at(multiSeq);
    std::vector<std::tuple<int, size_t, bool>> neighbors = mcp.getNeighbors();
    GridStorage& tStorage = grids.at(mcp.getDominateClass())->getStorage();
    std::cout << "neighbors: ";
    for (size_t n = 0 ; n < neighbors.size() ; n++) {
         std::cout << std::get<0>(neighbors.at(n)) << " ";
    }
    std::cout <<  std::endl;
    for (size_t n = 0 ; n < neighbors.size() ; n++) {
        int nextPoint = std::get<0>(neighbors.at(n));
        std::cout << "refine neighbor: " << nextPoint << std::endl;
        sgpp::base::MultipleClassPoint neighP = points.at(nextPoint);
        size_t dim = std::get<1>(neighbors.at(n));
        bool isLeft = std::get<2>(neighbors.at(n));
        
        GridStorage& nStorage = grids.at(neighP.getDominateClass())->getStorage();

        index_t source_index;
        level_t source_level;
        point.get(dim, source_level, source_index);

        // create the point
        // add points to classes if not yet inside the grids
        // adds original point to needed grid
        if ( !tStorage.isContaining(point) ) {
            // dominate class not set, insert point
            // should not happen
            std::cout << "!addP: " << multiSeq << " class: " << mcp.getDominateClass() << std::endl;
            createGridpoint(tStorage, point);
        }
        if ( !nStorage.isContaining(point) ) {
            // dominate class of neighbor not set, insert point
            std::cout << " addP: " << multiSeq << " class: " << neighP.getDominateClass() << std::endl;
            createGridpoint(nStorage, point);
        }
        // adds neighbor to needed grid
        GridPoint nPoint = multigrid.getStorage().getPoint(nextPoint);       
        if ( !tStorage.isContaining(nPoint) ) {
            // dominate class of neighbor not set, insert point
            createGridpoint(tStorage, nPoint);
            std::cout << " addP: " << nextPoint << " class: "
                    << mcp.getDominateClass() << std::endl;
        }
        if ( !nStorage.isContaining(nPoint) ) {
            // dominate class not set, insert point
            // should not happen
            createGridpoint(nStorage, nPoint);
            std::cout << "!addP: " << nextPoint << " class: "
                    << neighP.getDominateClass() << std::endl;
        }

        // create children in given dimension/direction
        // from refineGridpoint1D
        if ( isLeft ) {
            // generate left child, if necessary
            point.set(dim, source_level + 1, 2 * source_index - 1);
        } else {
            // generate right child, if necessary
            point.set(dim, source_level + 1, 2 * source_index + 1);
        }
        std::string pInfo = "(";
        for (size_t d = 0; d < multigrid.getStorage().getDimension(); d++) {
            pInfo += std::to_string(point.getStandardCoordinate(d)) + ",";
        }
        // insert point in both dominate classes
        if (!tStorage.isContaining(point)) {
            // point.setLeaf(true);
            // look for classes to insert
            std::cout << isLeft << " CaddP: " << pInfo << ") class: "
                    << mcp.getDominateClass() << " dim: " << dim << std::endl;
            createGridpoint(tStorage, point);
        }
        if (!nStorage.isContaining(point)) {
            // point.setLeaf(true);
            // look for classes to insert
            std::cout << isLeft << " CaddP: " << pInfo << ") class: "
                    << neighP.getDominateClass() << " dim: " << dim << std::endl;
            createGridpoint(nStorage, point);
        }
        point.set(dim, source_level, source_index);
    }
}

void HashRefinementMultipleClass::refineGridpointsCollection(GridStorage& storage,
            RefinementFunctor& functor,
            AbstractRefinement::refinement_container_type& collection) {
    // TODO (degel_kn): override source: HashRefinement
    // is this method needed?
    std::cout << " refineGridpointsCollection: " << std::endl;
    double threshold = functor.getRefinementThreshold();
    
    for (AbstractRefinement::refinement_pair_type& pair : collection) {
        if (pair.second >= threshold) {
            refineGridpoint(storage, pair.first->getSeq());
        }
    }
}

} /* namespace base */
} /* namespace sgpp */
