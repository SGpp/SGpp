// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "MultipleClassRefinementFunctor.hpp"


#include <iostream>
#include <tuple>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementMultipleClass.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>
#include <sgpp/base/tools/MultipleClassPoint.hpp>


namespace sgpp {
namespace datadriven {

MultipleClassRefinementFunctor::MultipleClassRefinementFunctor(std::vector<base::Grid*> grids,
                                std::vector<base::DataVector*> alphas,
                                size_t refinements_num, double thresh) :
                ZeroCrossingRefinementFunctor(grids, alphas, refinements_num,
                                false, false, thresh) {
    // create grid
    // prepareGrid();
    topPercent = 0.2;
    borderPenalty = 1;
}

double MultipleClassRefinementFunctor::operator()(base::GridStorage&
                              storage, size_t seq) const {
    // TODO (degel_kn): update scoring by class
    // score methode
    // check bool if use combined calc
    // find coresponing child in combined grid (class in grid_index) (class grid seq given)
    // ->getStorage().getSequenceNumber(gp)...
    // score for this point using combined grid infos
    // used only class related infos
    
    if (!refineMulti) {
        // use grid index for class, calc with combined grid
        base::HashGridPoint& gPoint = grids.at(current_grid_index)
                ->getStorage().getPoint(seq);
        size_t seqNew = multigrid->getStorage().getSequenceNumber(gPoint);
        seq = seqNew;
    }
    if ( seq >= points.size() ) {
        // new added Points not scored
        return 0.0;
    }
    
    base::MultipleClassPoint mcp = points.at(seq);
    base::HashGridPoint& gp = multigrid->getStorage().getPoint(seq);
    std::vector<std::tuple<size_t, size_t, bool>> neighbors = mcp.getNeighbors();

    double score = 0.0;
    // score neighbors
    for (size_t n = 0 ; n < neighbors.size() ; n++) {
        size_t nextPoint = std::get<0>(neighbors.at(n));
        size_t dimension = std::get<1>(neighbors.at(n));

        // change in class
        double ctdt = mcp.getDensity(mcp.getDominateClass());
        double ctdn = mcp.getDensity(points.at(nextPoint).getDominateClass());
        double cndt = points.at(nextPoint).getDensity(mcp.getDominateClass());
        double cndn = points.at(nextPoint).getDensity(points.at(nextPoint).getDominateClass());

        double sCon = std::abs((ctdt - cndt) - (ctdn - cndn)) /
                static_cast<double>(grids.at(current_grid_index)->getDimension());

        // level penalising
        size_t maxD = multigrid->getStorage().getPoint(nextPoint).getLevel(dimension);
        maxD = ( maxD < gp.getLevel(dimension) ) ? gp.getLevel(dimension) : maxD;

        sCon = sCon / pow(2.0, maxD);
        if (refineMulti) {
            score = score + sCon;
        } else {
            bool c1 = mcp.getDominateClass() == current_grid_index;
            bool c2 = points.at(nextPoint).getDominateClass() == current_grid_index;
            if ( c1 || c2 ) {
                // one dominate class is scored
                score = score + sCon;
            }
        }
    }

    // score borders
    double scoreB = borderPenalty * (mcp.getBorderScore() /
                    static_cast<double>(grids.at(0)->getDimension()));
    if (refineMulti) {
        scoreB = scoreB * mcp.getDensity(mcp.getDominateClass());
    } else {
        scoreB = scoreB * mcp.getDensity(current_grid_index);
    }
    if (scoreB > 0) {
        score = score + scoreB;
    }
    borderCnt += 1.0;
    borderSum += mcp.getBorderScore() * mcp.getDensity(mcp.getDominateClass());

    // score point on its own
    std::vector<std::tuple<double, size_t, bool>> top = mcp.getTopClasses(topPercent);
    if (refineMulti) {
        score = score * (1+(top.size()/static_cast<double>(grids.size())));
    } else {
        for (size_t n = 0 ; n < top.size() ; n++) {
            if (std::get<1>(top.at(n)) == current_grid_index) {
                // add to score if class in close classes
                score = score * (1+(top.size()/static_cast<double>(grids.size())));
                break;
            }
        }
    }
    return score;
}

void MultipleClassRefinementFunctor::prepareGrid() {
    size_t dim = grids.at(0)->getDimension();
    points.clear();

    multigrid = base::Grid::createLinearGrid(dim);
    for (size_t n = 0; n < grids.size() ; n++) {
        for (size_t seq = 0; seq < grids.at(n)->getSize(); seq++) {
            base::HashGridPoint& gp = grids.at(n)->getStorage().getPoint(seq);
            if ( !multigrid->getStorage().isContaining(gp) ) {
                // insert a new point
                unsigned int* level = new unsigned int[dim];
                unsigned int* index = new unsigned int[dim];
                bool leaf = gp.isLeaf();
                for (size_t d = 0; d < dim; d++) {
                    level[d] = gp.getLevel(d);
                    index[d] = gp.getIndex(d);
                }
                multigrid->insertPoint(dim, level, index, leaf);

                points.push_back(base::MultipleClassPoint(gp,
                        grids, alphas));
            }
        }
    }
    multigrid->getStorage().recalcLeafProperty();
    findCrossings(-1, -1, 0, 0);
}

void MultipleClassRefinementFunctor::findCrossings(
            size_t leftP, size_t rightP, size_t seq, size_t dim) {
    if ( seq >= points.size() ) {
        return;
    }
    base::HashGridPoint& gp = multigrid->getStorage().getPoint(seq);

    // Find the geometric neighbors
    base::HashGridIterator iter(multigrid->getStorage());
    for (size_t d = 0; d < multigrid->getStorage().getDimension(); d++) {
      // Left neighbor
      if ( hasChild(gp, d, true) ) {
        base::HashGridPoint child = base::HashGridPoint(gp);
        getChild(gp, d, true, child);
        if ( d != dim ) {
            findCrossings(-1, seq,
                        multigrid->getStorage().getSequenceNumber(child), d);
        } else {
            findCrossings(leftP, seq,
                        multigrid->getStorage().getSequenceNumber(child), d);
        }
      } else {
        if ( d == dim && seq < points.size() && leftP < points.size() ) {
                // std::cout << " neighbor: "<< leftP << std::endl;
            if (points.at(seq).getDominateClass() !=
                        points.at(leftP).getDominateClass()) {
                points.at(leftP).addNeighbor(seq, d, false);
                points.at(seq).addNeighbor(leftP, d, true);
            }
        }
        if ( gp.getIndex(d) > 1 ) {
        } else {
            // add border
            points.at(seq).addBorder(d, gp.getLevel(d), true);
        }
      }

      // Right neighbor
      if ( hasChild(gp, d, false) ) {
        base::HashGridPoint child = base::HashGridPoint(gp);
        getChild(gp, d, false, child);
        if ( d != dim ) {
            findCrossings(seq, -1,
                    multigrid->getStorage().getSequenceNumber(child), d);
        } else {
            findCrossings(seq, rightP,
                    multigrid->getStorage().getSequenceNumber(child), d);
        }
      } else {
        // std::cout << seq << ": l neighbor in dim: " << d;
        if ( d == dim && seq < points.size() && rightP < points.size() ) {
            // has neighbor in direction
            // std::cout << " neighbor: "<< rightP << std::endl;
            if ( points.at(seq).getDominateClass() !=
                        points.at(rightP).getDominateClass() ) {
                points.at(seq).addNeighbor(rightP, d, false);
                points.at(rightP).addNeighbor(seq, d, true);
            }
        }
        if ( gp.getIndex(d) < pow(2.0, gp.getLevel(d)) - 1 ) {
        } else {
            // add border
            points.at(seq).addBorder(d, gp.getLevel(d), false);
        }
      }
    }
}

void MultipleClassRefinementFunctor::refine(size_t partCombined) {
    prepareGrid();

    base::HashRefinementMultipleClass refine(*multigrid, &points, grids,
            borderSum, borderCnt, topPercent);
    borderSum = 0.0;
    borderCnt = 0.0;

    // refine every class for itself
    size_t tmp = refinements_num;
    refineMulti = false; // use grid index
    refinements_num = refinements_num - partCombined;
    if (refinements_num > 0) {
        for (size_t i = 0; i < grids.size() ; i++) {
            setGridIndex(i);
            refine.free_refine(grids.at(i)->getStorage(), *this);
            borderSum = 0.0;
            borderCnt = 0.0;
        }
    }

    // refine combined grid
    refineMulti = true;
    refinements_num = partCombined;
    if (refinements_num > 0) {
        refine.free_refine(multigrid->getStorage(), *this);
    }

    // printScores();
    // printPointsPlott();
    // printPointsInfo();
    refinements_num = tmp;
}

double MultipleClassRefinementFunctor::getTopPercent() {
    return topPercent;
}

void MultipleClassRefinementFunctor::setTopPercent(double newPercent) {
    topPercent = newPercent;
}

double MultipleClassRefinementFunctor::getBorderPenalty() {
    return borderPenalty;
}

void MultipleClassRefinementFunctor::setBorderPenalty(double newPenalty) {
    borderPenalty = newPenalty;
}

bool MultipleClassRefinementFunctor::hasChild(const base::HashGridPoint& gp,
              size_t d, bool left) const {
    base::HashGridIterator iter(multigrid->getStorage());
    iter.set(gp);
    if ( left ) {
      return iter.hintLeft(d);
    } else {
      return iter.hintRight(d);
    }
}

void MultipleClassRefinementFunctor::printPointsPlott() {
    prepareGrid();
    for (size_t i = 0; i < points.size() ; i++) {
        // used as name for printing
        std::cout << i << ",";
        // coordinates in all dimensions
        sgpp::base::GridPoint& gp = multigrid->getStorage().getPoint(i);
        for (size_t d = 0; d < multigrid->getStorage().getDimension(); d++) {
            std::cout << gp.getStandardCoordinate(d) << ",";
        }
        sgpp::base::MultipleClassPoint mcp = points.at(i);
        for (size_t d = 0; d < grids.size(); d++) {
            std::cout << mcp.isClassSet(d) << ",";
        }
        // dominate class
        std::cout << points.at(i).getDominateClass() << std::endl;
    }
}

void MultipleClassRefinementFunctor::printPointsInfo() {
    prepareGrid();
    for (size_t i = 0; i < points.size() ; i++) {
        std::cout << "Point " << i << ": ";
        std::cout << points.at(i).getDominateClass() << " ~";

        sgpp::base::GridPoint& gp = multigrid->getStorage().getPoint(i);
        for (size_t d = 0; d < multigrid->getStorage().getDimension(); d++) {
            std::cout << " " << gp.getStandardCoordinate(d) << " -";
        }
        std::cout << points.at(i).toString() << std::endl;
    }
}
} /* namespace datadriven */
} /* namespace sgpp */
