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
#include <sgpp/datadriven/application/MultipleClassPoint.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementMultipleClass.hpp>


namespace sgpp {
namespace datadriven {

MultipleClassRefinementFunctor::MultipleClassRefinementFunctor(std::vector<base::Grid*> grids,
                                std::vector<base::DataVector*> alphas,
                                size_t refinements_num, bool level_penalize,
                                bool pre_compute, double thresh) :
    ZeroCrossingRefinementFunctor(grids, alphas, refinements_num,
                    level_penalize, pre_compute, thresh){
    // create grid
    prepareGrid(grids, alphas);
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
        size_t seqNew = grids.at(current_grid_index)
                ->getStorage().getSequenceNumber(gPoint);
        seq = seqNew;
    }
    
    
    sgpp::datadriven::MultipleClassPoint mcp = points.at(seq);
    base::HashGridPoint& gp = multigrid->getStorage().getPoint(seq);
    std::vector<std::tuple<int, size_t, bool>> neighbors = mcp.getNeighbors();

    double score = 0.0;
    /* redesign for class dependent score
    // score point on its own
    std::vector<std::tuple<double, int, bool>> top = mcp.getTopClasses(0.2);
    double score = (static_cast<double>(top.size())-1.0) * 0.05 / static_cast<double>(grids.size());
    // std::cout  << seq << " classes top: " << score << std::endl; */

    // score connections
    for (size_t n = 0 ; n < neighbors.size() ; n++) {
        int nextPoint = std::get<0>(neighbors.at(n));
        size_t dimension = std::get<1>(neighbors.at(n));

        // score on change in class
        // std::cout << "class change: " << std::endl;
        double dom1a = mcp.getDensity(mcp.getDominateClass());
        // std::cout << "     dom1a: " << dom1a << std::endl;
        double dom2a = mcp.getDensity(points.at(nextPoint).getDominateClass());
        // std::cout << "     dom2a: " << dom2a << std::endl;
        double dom1b = points.at(nextPoint).getDensity(mcp.getDominateClass());
        // std::cout << "     dom1b: " << dom1b << std::endl;
        double dom2b = points.at(nextPoint).getDensity(points.at(nextPoint).getDominateClass());
        // std::cout << "     dom2b: " << dom2b << std::endl;
        // double test1 = (dom1a - dom1b) + (dom2a - dom2b);
        // std::cout << "std::abs test1: " << std::abs(test1) << std::endl;
        // test2 for class dependent scoring v1
        double test2 = (dom1a - dom1b) - (dom2a - dom2b);
        // std::cout << "test2: " << test2 << std::endl;
        double test3 = (dom1a + dom2b) + (dom1a - dom1b) - (dom2a - dom2b);
        // double test4 = (dom1a) * (dom1a - dom1b) - (dom2b) * (dom2a - dom2b);
        // std::cout << "test4: " << test4 << std::endl;
        double sCon = test2 / static_cast<double>(grids.at(current_grid_index)->getDimension());

        // score on level
        size_t maxD = multigrid->getStorage().getPoint(nextPoint).getLevel(dimension);
        // std::cout << "neighbors level: " << gp.getLevel(dimension) << "~"<< maxD;
        maxD = ( maxD < gp.getLevel(dimension) ) ? gp.getLevel(dimension) : maxD;
        // std::cout << ": max  " << (1.0 /maxD) << std::endl;
        sCon = sCon / static_cast<double>(maxD);

        if (!refineMulti) {
            bool c1 = mcp.getDominateClass() == current_grid_index;
            bool c2 = points.at(nextPoint).getDominateClass() == current_grid_index;
            if ( c1 || c2 ) {
                // one dominate class is looked for
                score = score + sCon;
            }
        } else {
            score = score + sCon;
        }
        // std::cout << "con " << seq<< "-"<< nextPoint <<": ";
        // std::cout << sCon << "    levelDim  " << (1.0 /maxD) << std::endl;
    }
    if (refineMulti) {
        scoresToPrint.at(seq).append("m");
    } else {
        scoresToPrint.at(seq).append(std::to_string(current_grid_index));
    }
    scoresToPrint.at(seq).append(": " + std::to_string(score) + " | ");

    // std::cout << "Point: " << seq << " end score: " << score <<  std::endl;
    return score;
}

void MultipleClassRefinementFunctor::prepareGrid(
                                std::vector<base::Grid*> gridsNew,
                                std::vector<base::DataVector*> alphasNew) {
    grids = gridsNew;
    alphas = alphasNew;
    size_t dim = grids.at(0)->getDimension();
    points.clear();

    multigrid = sgpp::base::Grid::createLinearGrid(dim);
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

                points.push_back(sgpp::datadriven::MultipleClassPoint(gp,
                        grids, alphas));
            }
        }
    }
    multigrid->getStorage().recalcLeafProperty();
    size_t allP = multigrid->getStorage().getSize() + 10;
    findCrossings(allP, allP, 0, 0);

    if ( scoresToPrint.size() != points.size() ) {
        scoresToPrint.clear();
        std::vector<std::string> v(points.size(),  std::string(""));
        scoresToPrint = v;
    }
}

void MultipleClassRefinementFunctor::findCrossings(
            size_t leftP, size_t rightP, size_t seq, size_t dim) {

    if ( seq > points.size() ) {
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
            if (points.at(seq).getDominateClass() !=
                        points.at(leftP).getDominateClass()) {
                points.at(leftP).addNeighbor(seq, d, true);
                points.at(seq).addNeighbor(leftP, d, false);
            }
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
        if ( d == dim && seq < points.size() && rightP < points.size() ) {
            if ( points.at(seq).getDominateClass() !=
                        points.at(rightP).getDominateClass() ) {
                points.at(seq).addNeighbor(rightP, d, true);
                points.at(rightP).addNeighbor(seq, d, false);
            }
        }
      }
    }
}

void MultipleClassRefinementFunctor::refine() {
    // Settings: local for now
    // points to refine in combine grid, smaller than refinements_num
    int partComined = 0; // 0-3 for now
    // TODO (degel_kn): update refine generall?
    // refine methode
    // insert in both combined and classgrid
    // check if neighbors have given class, insert neighbors

    // code from generator.refine():
    /*  HashRefinement refine;
        refine.free_refine(this->storage, func);*/
    // update for HashRefinementMultipleClass
    base::HashRefinementMultipleClass refine;

    // refine every class for itself
    int tmp = refinements_num;
    refineMulti = false; // use grid index
    std::cout << "num ref: " << refinements_num << std::endl;
    refinements_num = refinements_num - partComined;
    std::cout << "num comb: " << refinements_num << std::endl;
    for (size_t i = 0; i < grids.size() ; i++) {
        setGridIndex(i);
        // grids.at(i)->getGenerator().refine(*this);
        refine.free_refine(grids.at(i)->getStorage(), *this);
    }

    // refine combined grid
    refineMulti = true;
    refinements_num = partComined;
    refine.free_refine(multigrid->getStorage(), *this);

    printPointsInfo();
    printScores();
    refinements_num = tmp;
}

void MultipleClassRefinementFunctor::refineCombinedGrid() {
    // TODO (degel_kn): seperate mulitgrid?
}

base::Grid* MultipleClassRefinementFunctor::getCombinedGrid() {
    return multigrid;
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
    for (size_t i = 0; i < points.size() ; i++) {
        // used as name for printing
        std::cout << i << ",";
        // coordinates in all dimensions
        sgpp::base::GridPoint& gp = multigrid->getStorage().getPoint(i);
        for (size_t d = 0; d < multigrid->getStorage().getDimension(); d++) {
            std::cout << gp.getStandardCoordinate(d) << ",";
        }
        // dominate class
        std::cout << points.at(i).getDominateClass() << std::endl;
    }
}

void MultipleClassRefinementFunctor::printPointsInfo() {
    for (size_t i = 0; i < points.size() ; i++) {
        std::cout << "Point " << i << ": ";
        std::cout << points.at(i).getDominateClass() << " ~ ";

        sgpp::base::GridPoint& gp = multigrid->getStorage().getPoint(i);
        for (size_t d = 0; d < multigrid->getStorage().getDimension(); d++) {
            std::cout << gp.getStandardCoordinate(d) << " - ";
        }
        std::cout << points.at(i).toString() << std::endl;
    }
}

void MultipleClassRefinementFunctor::printScores() {
    for (size_t i = 0; i < scoresToPrint.size() ; i++) {
        std::cout << i << " - "<< scoresToPrint.at(i) << std::endl;
    }
    
}

} /* namespace datadriven */
} /* namespace sgpp */
