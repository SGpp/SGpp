// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "MultipleClassRefinementFunctor.hpp"

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/application/MultipleClassPoint.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <tuple>
#include <cmath>
#include <unordered_set>


namespace sgpp {
namespace datadriven {

MultipleClassRefinementFunctor::MultipleClassRefinementFunctor(std::vector<base::Grid*> grids,
                                std::vector<base::DataVector*> alphas,
                                size_t refinements_num,
                                bool level_penalize,
                                bool pre_compute,
                                double thresh) :
    ZeroCrossingRefinementFunctor(grids,
                                alphas,
                                refinements_num,
                                level_penalize,
                                pre_compute,
                                thresh){
    base::DataVector coords(grids.at(0)->getDimension());
        
    //create grid
    std::vector<std::unique_ptr<sgpp::base::OperationEval>> evalOps;
    size_t dim = grids.at(0)->getDimension();
    sgpp::base::DataVector p(dim);
    for (size_t i = 0; i < grids.size(); i++) {
        std::unique_ptr<sgpp::base::OperationEval>
            e(sgpp::op_factory::createOperationEval(*grids.at(i)));
        evalOps.push_back(std::move(e));
    }
    double v = 0.0;

    multigrid = sgpp::base::Grid::createLinearGrid(dim);
    for (int n = 0; n < grids.size() ; n++) {
        for (int seq = 0; seq < grids.at(n)->getSize(); seq++) {
            base::HashGridPoint& gp = grids.at(n)->getStorage().getPoint(seq);
            if (! multigrid->getStorage().isContaining(gp)) {
                // insert a new point
                unsigned int* level = new unsigned int[dim];
                unsigned int* index = new unsigned int[dim];
                bool leaf = gp.isLeaf();
                for (size_t d = 0; d < dim; d++) {
                    level[d] = gp.getLevel(d);
                    index[d] = gp.getIndex(d);
                }
                multigrid->insertPoint(dim, level, index, leaf);

                points.push_back(sgpp::datadriven::MultipleClassPoint(gp, grids, alphas));
                /*for (int t =  0 ; t < grids.size() ; t++) {
                    std::vector<double> vec = getEvalVector(t ,seq); 
                    v = evalOps.at(t)->eval(*alphas.at(t), p);
                    //v = vec.at(t);
                    // TODO calc useful
                    gp.getStandardCoordinates(coords);
                    points.back().updateClass( t , v , grids.at(t)->getStorage().isContaining(gp) );
                }
                points.back().resortClasses();*/
            }
            
        }
    }
    multigrid->getStorage().recalcLeafProperty();
    
    findCrossings(-1, -1, 0, 0);
    
    std::cout << "Points: " << multigrid->getStorage().getSize() << std::endl;
    std::cout << "Points: " << points.size() << std::endl;
    for (size_t i = 0; i < points.size(); i++) {
        sgpp::base::GridPoint& gp = multigrid->getStorage().getPoint(i);

        std::cout << "Point " << i << ": " << points.at(i).getDominateClass() << " leaf " << gp.isLeaf();
        std::cout << " ~ " << gp.getStandardCoordinate(0) << " - " << gp.getStandardCoordinate(1);
        std::cout << std::endl <<  points.at(i).toString() << std::endl;
    }
    std::cout << std::endl << std::endl;
    std::cout << "Print for plotting:"<< std::endl;
    for (size_t i = 0; i < points.size() ; i++) {
        sgpp::base::GridPoint& gp = multigrid->getStorage().getPoint(i);
        std::cout << i << ",";
        std::cout << gp.getStandardCoordinate(0) << "," << gp.getStandardCoordinate(1);
        std::cout << "," << points.at(i).getDominateClass() << std::endl;
    }
    std::cout << std::endl << std::endl;


}

double MultipleClassRefinementFunctor::operator()(base::GridStorage&
                                                   storage,
                                                   size_t seq) const {
    sgpp::datadriven::MultipleClassPoint mcp = points.at(seq);
    base::HashGridPoint& gp = multigrid->getStorage().getPoint(seq);
    std::vector<std::tuple<int, int>> neighbors = mcp.getNeighbors();

    // score point on its own
    std::vector<std::tuple<double, int, bool>> top = mcp.getTopClasses(0.2);
    double score = (top.size()-1) * 0.05 / grids.size();
    std::cout  << seq << " classes top: " << score << std::endl;

    // score connections
    for (int n = 0 ; n < neighbors.size() ; n++) {
        int nextPoint = std::get<0>(neighbors.at(n));
        int dimension = std::get<1>(neighbors.at(n));

        //score on change in class
        //std::cout << "class change: " << std::endl;
        double dom1a = mcp.getDensity(mcp.getDominateClass());
        //std::cout << "     dom1a: " << dom1a << std::endl;
        double dom2a = mcp.getDensity(points.at(nextPoint).getDominateClass());
        //std::cout << "     dom2a: " << dom2a << std::endl;
        double dom1b = points.at(nextPoint).getDensity(mcp.getDominateClass());
        //std::cout << "     dom1b: " << dom1b << std::endl;
        double dom2b = points.at(nextPoint).getDensity(points.at(nextPoint).getDominateClass());
        //std::cout << "     dom2b: " << dom2b << std::endl;
        double test1 = (dom1a - dom1b) + (dom2a - dom2b);
        //std::cout << "std::abs test1: " << std::abs(test1) << std::endl;
        double test2 = (dom1a - dom1b) - (dom2a - dom2b);
        //std::cout << "test2: " << test2 << std::endl;
        double test3 = (dom1a + (dom1a - dom1b)) - (dom2b +(dom2a - dom2b));
        //std::cout << "test3: " << test3 << std::endl;
        double test4 = (dom1a) * (dom1a - dom1b) - (dom2b) * (dom2a - dom2b);
        //std::cout << "test4: " << test4 << std::endl;
        
        double sCon = test4 / grids.at(0)->getDimension();
        
        //score on level
        int maxD = multigrid->getStorage().getPoint(nextPoint).getLevel(dimension);
        maxD = (maxD>gp.getLevel(dimension))? gp.getLevel(dimension) : maxD ;
        //std::cout << "neighbors level: max  " << (1.0 /maxD) << " dim: " << maxD << std::endl;
        sCon = sCon / maxD;
        
        score = score + sCon;
        std::cout << "con " << seq<< "-"<< nextPoint <<": "<< sCon << "    levelDim  " << (1.0 /maxD) << std::endl;
    }

    std::cout << "Point: " << seq << " end score: " << score <<  std::endl;
    return score;
}

void MultipleClassRefinementFunctor::findCrossings(
            int leftP, int rightP, int seq, size_t dim) {


    if (seq < 0  || seq > points.size() ) {
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
            findCrossings(-1, seq, multigrid->getStorage().getSequenceNumber(child), d);
        }else {
            findCrossings(leftP, seq, multigrid->getStorage().getSequenceNumber(child), d);
        }
      } else {
        if ( d == dim ) {
            if ( seq >= 0  || seq < points.size() ) {
                if ( leftP >= 0  || leftP < points.size() ) {
                    if (points.at(seq).getDominateClass() != points.at(leftP).getDominateClass()) {
                        std::cout << d << "-Class change: " << seq << " <-> " << leftP << std::endl;
                        points.at(leftP).addNeighbor(seq, d);
                        points.at(seq).addNeighbor(leftP, d);
                    }
                }
            }
         }
      }

      // Right neighbor
      if (hasChild(gp, d, false)) {
        base::HashGridPoint child = base::HashGridPoint(gp);
        getChild(gp, d, false, child);
        if (d != dim) {
            findCrossings(seq, -1, multigrid->getStorage().getSequenceNumber(child), d);
        } else {
            findCrossings(seq, rightP, multigrid->getStorage().getSequenceNumber(child), d);
        }
      } else {
        if (d == dim) {
            if ( seq >= 0  || seq < points.size() ) {
                if ( rightP >= 0  || rightP < points.size() ) {
                    if (points.at(seq).getDominateClass() != points.at(rightP).getDominateClass()) {
                        std::cout << d << "-Class change: " << seq << " <-> " << rightP << std::endl;
                        points.at(seq).addNeighbor(rightP, d);
                        points.at(rightP).addNeighbor(seq, d);
                    }
                }
            }
         }
      }
    }
}

base::Grid* MultipleClassRefinementFunctor::getCombinedGrid(){
    return multigrid;
}

bool MultipleClassRefinementFunctor::hasChild(const base::HashGridPoint& gp,
                                               size_t d,
                                               bool left) const {
    base::HashGridIterator iter(multigrid->getStorage());
    iter.set(gp);
    if (left) {
      return iter.hintLeft(d);
    } else {
      return iter.hintRight(d);
    }
  }

} /* namespace datadriven */
} /* namespace sgpp */
