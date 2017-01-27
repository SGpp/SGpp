// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "MultipleClassRefinementFunctor.hpp"

namespace sgpp {
namespace datadriven {

MultipleClassRefinementFunctor::MultipleClassRefinementFunctor(std::vector<base::Grid*> grids,
                                std::vector<base::DataVector*> alphas,
                                std::vector<sgpp::datadriven::MultipleClassPoint> * pts,
                                base::GridStorage& store,
                                size_t refinements_num,
                                bool level_penalize,
                                bool pre_compute,
                                double thresh) :
    ZeroCrossingRefinementFunctor(grids,
                                alphas,
                                refinements_num,
                                level_penalize,
                                pre_compute,
                                thresh), points(pts), storage(store) {

    findCrossings(-1, -1, 0, 0);
}

double MultipleClassRefinementFunctor::operator()(base::GridStorage&
                                                   storage,
                                                   size_t seq) const {
    sgpp::datadriven::MultipleClassPoint mcp = points->at(seq);
    std::vector<int> neighbors = mcp.getNeighbors();

    // score point on its own
    std::vector<std::tuple<double, int, bool>> top = mcp.getTopClasses(0.1);
    double score = (top.size()-1) * 0.5;

    // score connections
    for (int n = 0; n < neighbors.size() ; n++) {
        if (neighbors.at(n) < points->size() && neighbors.at(n) >= 0) {
            if (mcp.getDominateClass() != points->at(neighbors.at(n)).getDominateClass()){
                score = score + 1;
            }
        } 
    }

    std::cout << "Point: " << seq << " score: " << score << std::endl;
    return score;
}

void MultipleClassRefinementFunctor::findCrossings(
            int leftP, int rightP, int seq, size_t dim) {


    if (seq < 0  || seq > points->size() ) {
        return;
    }
    base::HashGridPoint& gp = storage.getPoint(seq);

    // Find the geometric neighbors
    base::HashGridIterator iter(storage);
    for (size_t d = 0; d < storage.getDimension(); d++) {
      // Left neighbor
      if ( hasChild(gp, d, true) ) {
        base::HashGridPoint child = base::HashGridPoint(gp);
        getChild(gp, d, true, child);
        if ( d != dim ) {
            findCrossings(-1, seq, storage.getSequenceNumber(child), d);
        }else {
            findCrossings(leftP, seq, storage.getSequenceNumber(child), d);
        }
      } else {
          //addNeigbor(leftP, seq, d);
        if ( d == dim ) {
            if ( seq >= 0  || seq < points->size() ) {
                points->at(seq).addNeighbor(leftP);
                if ( leftP >= 0  || leftP < points->size() ) {
                    if (points->at(seq).getDominateClass() != points->at(leftP).getDominateClass()) {
                        std::cout << "Class change: " << seq << " <-> " << leftP << std::endl;
                    }
                }
            }
            if ( leftP >= 0  || leftP < points->size() ) {
                points->at(leftP).addNeighbor(seq);
            }
         }
      }

      // Right neighbor
      if (hasChild(gp, d, false)) {
        base::HashGridPoint child = base::HashGridPoint(gp);
        getChild(gp, d, false, child);
        if (d != dim) {
            findCrossings(seq, -1, storage.getSequenceNumber(child), d);
        } else {
            findCrossings(seq, rightP, storage.getSequenceNumber(child), d);
        }
      } else {
          //addNeigbor(seq, rightP, d);
        if (d == dim) {
            if ( seq >= 0  || seq < points->size() ) {
                points->at(seq).addNeighbor(rightP);
                if ( rightP >= 0  || rightP < points->size() ) {
                    if (points->at(seq).getDominateClass() != points->at(rightP).getDominateClass()) {
                        std::cout << "Class change: " << seq << " <-> " << rightP << std::endl;
                    }
                }
            }
             
            if ( rightP >= 0  || rightP < points->size() ) {
                points->at(rightP).addNeighbor(seq);
            }
         }
      }
    }
    /*
    void MultipleClassRefinementFunctor::addNeigbor(int leftP, int rightP, int dim) {
        sgpp::datadriven::MultipleClassPoint* left = nullptr;
        sgpp::datadriven::MultipleClassPoint* right = nullptr;
        if (leftP >= 0  || leftP < points->size() ) {
            left = points->at(leftP);
        }
        if (rightP >= 0  || rightP < points->size() ) {
            right = points->at(rightP)
        }
        if (left) {
            left->addNeighbor(right, dim, false);
        }
        if (right) {
            right->addNeighbor(left, dim, true);
        }
    }*/
}

} /* namespace datadriven */
} /* namespace sgpp */
