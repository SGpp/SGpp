// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <vector>
#include <map>
#include <string>


namespace sgpp {
namespace datadriven {

  ZeroCrossingRefinementFunctor::
  ZeroCrossingRefinementFunctor(std::vector<base::Grid*> grids,
                                std::vector<base::DataVector*> alphas,
                                std::vector<double> priors,
                                size_t refinements_num,
                                bool level_penalize,
                                bool pre_compute,
                                double thresh) :
    grids(grids), alphas(alphas), priors(priors), current_grid_index(0),
    refinements_num(refinements_num), threshold(thresh),
    level_penalize(level_penalize),
    pre_compute(pre_compute), pre_comp_evals() {
    for (size_t i = 0; i < grids.size(); i++) {
      this->pre_comp_evals.push_back(std::map<std::string, double>());
    }
  }

  double ZeroCrossingRefinementFunctor::operator()(base::GridStorage&
                                                   storage,
                                                   size_t seq) const {
    base::HashGridPoint& gp = storage.getPoint(seq);
    std::vector<double> seqEvals = getEvalVector(current_grid_index, seq);
    double levelSum = gp.getLevelSum();
    double levelW = pow(2.0, -levelSum);

    // Find the geometric neighbors
    base::HashGridIterator iter(storage);
    std::vector<size_t> neighSeq;
    for (size_t d = 0; d < storage.getDimension(); d++) {
      // Left neighbor
      if (hasChild(gp, d, true)) {
        base::HashGridPoint child = base::HashGridPoint(gp);
        base::HashGridPoint down = base::HashGridPoint(gp);
        getChild(gp, d, true, child);
        goDown(child, down, d, false);
        neighSeq.push_back(storage.getSequenceNumber(down));
      } else {
        // Check if left-most grid-point on level, which has no left neigh
        if (gp.getIndex(d) > 1) {
          base::HashGridPoint gp_c = base::HashGridPoint(gp);
          base::HashGridPoint up = base::HashGridPoint(gp);
          goUp(gp_c, up, d, true);
          neighSeq.push_back(storage.getSequenceNumber(up));
        }
      }
      // Right neighbor
      if (hasChild(gp, d, false)) {
        base::HashGridPoint child = base::HashGridPoint(gp);
        base::HashGridPoint down = base::HashGridPoint(gp);
        getChild(gp, d, false, child);
        goDown(child, down, d, true);
        neighSeq.push_back(storage.getSequenceNumber(down));
      } else {
        // Check if right-most grid point on level, which has no right neigh
        if (gp.getIndex(d) < pow(2.0, gp.getLevel(d)) - 1) {
          base::HashGridPoint gp_c = base::HashGridPoint(gp);
          base::HashGridPoint up = base::HashGridPoint(gp);
          goUp(gp_c, up, d, false);
          neighSeq.push_back(storage.getSequenceNumber(up));
        }
      }
    }

    std::vector<double> neighEvals;
    double maxScore = 0.0;
    // Iterate over all (gridpoint, neighbor) pairs and compute the score
    for (size_t i = 0; i < neighSeq.size(); i++) {
      neighEvals = getEvalVector(current_grid_index, neighSeq.at(i));
      // eval at seq, current grid
      double currThisEval = seqEvals.at(current_grid_index);
      // eval at neigh, current grid
      double currNeighEval = neighEvals.at(current_grid_index);
      for (size_t j = 0; j < grids.size(); j++) {
        if (j == current_grid_index) {
          continue;
        }
        // eval at seq, j-th grid
        double thisEval = seqEvals.at(j);
        // eval at neigh, j-th grid
        double neighEval = neighEvals.at(j);
        if (sgn(currThisEval - thisEval) != sgn(currNeighEval - neighEval)) {
          double diff = fabs(currThisEval - thisEval);
          double neighDiff = fabs(currNeighEval - neighEval);
          double score = sqrt(diff * neighDiff) *
            fabs((currThisEval * thisEval));
          if (maxScore < score) {
            maxScore = score;
          }
        }
      }
    }
    if (level_penalize) {
      maxScore *= levelW;
    }
    return maxScore;
  }

  // Gets evaluations of grid point seq (at grid ind) at all grids
  // e.g. grid point with seq 4 at grid=0 is evaluated at all grids
  std::vector<double> ZeroCrossingRefinementFunctor::
  getEvalVector(size_t ind,
                size_t seq) const {
    base::HashGridPoint& gp = grids.at(ind)->getStorage().getPoint(seq);
    base::DataVector coords(grids.at(ind)->getDimension());
    gp.getStandardCoordinates(coords);
    std::vector<double> evals;
    if (pre_compute) {
      for (size_t i = 0; i < grids.size(); i++) {
        std::string key = coords.toString();
        evals.push_back(pre_comp_evals.at(i).at(key));
      }
    } else {
      for (size_t j = 0; j < grids.size(); j++) {
        std::unique_ptr<base::OperationEval>
          opEval(op_factory::createOperationEval(*grids.at(j)));
        double eval = (opEval->eval(*alphas.at(j), coords))*priors.at(j);
        evals.push_back(eval);
      }
    }
    return evals;
  }

  double ZeroCrossingRefinementFunctor::start() const {
    return 0.0;
  }

  size_t ZeroCrossingRefinementFunctor::getRefinementsNum() const {
    return this->refinements_num;
  }

  double ZeroCrossingRefinementFunctor::getRefinementThreshold() const {
    return this->threshold;
  }

  void ZeroCrossingRefinementFunctor::setGridIndex(size_t grid_index) {
    this->current_grid_index = grid_index;
  }

  size_t ZeroCrossingRefinementFunctor::getNumGrids() {
    return this->grids.size();
  }


  // For comments see GridPointBasedRefinementFunctor.cpp, exactly the same
  void ZeroCrossingRefinementFunctor::preComputeEvaluations() {
    base::DataVector p(grids.at(0)->getDimension());
    std::string key = "";
    double v = 0.0;
    for (size_t i = 0; i < grids.size(); i++) {
      std::unique_ptr<base::OperationEval>
        opEval(op_factory::createOperationEval(*grids.at(i)));
      pre_comp_evals.at(i).clear();

      for (size_t j = 0; j < grids.size(); j++) {
        for (size_t k = 0; k < grids.at(j)->getSize(); k++) {
          grids.at(j)->getStorage().getPoint(k).getStandardCoordinates(p);
          key = p.toString();
          if (pre_comp_evals.at(i).count(key) == 0) {
            v = (opEval->eval(*alphas.at(i), p))*priors.at(i);
            pre_comp_evals.at(i).insert(std::pair<std::string, double>(key,
                                                                       v));
          }
        }
      }
    }
  }

  int ZeroCrossingRefinementFunctor::sgn(double d) const {
    return d < 0 ? -1 : 1;
  }

  // Used to get geom. neighbor of non-leaf grid points
  void ZeroCrossingRefinementFunctor::goDown(base::HashGridPoint& gp,
                                             base::HashGridPoint& down,
                                             size_t d,
                                             bool left) const {
    // Child in direction exists? If not stop.
    if (hasChild(gp, d, left)) {
      getChild(gp, d, left, gp);
      goDown(gp, down, d, left);
    } else {
      unsigned int l = gp.getLevel(d);
      unsigned int i = gp.getIndex(d);
      down.set(d, l, i);
    }
  }

  // Use to get geom. neightbor of leaf grid point
  void ZeroCrossingRefinementFunctor::goUp(base::HashGridPoint& gp,
                                           base::HashGridPoint& up,
                                           size_t d,
                                           bool left) const {
    if (isLeftChild(gp, d) != left) {
      getParent(gp, d, up);
    } else {
      getParent(gp, d, gp);
      goUp(gp, up, d, left);
    }
  }

  bool ZeroCrossingRefinementFunctor::hasChild(const base::HashGridPoint& gp,
                                               size_t d,
                                               bool left) const {
    base::HashGridIterator iter(grids.at(current_grid_index)->getStorage());
    iter.set(gp);
    if (left) {
      return iter.hintLeft(d);
    } else {
      return iter.hintRight(d);
    }
  }

  bool ZeroCrossingRefinementFunctor::isLeftChild(const base::HashGridPoint& gp,
                                                  size_t d) const {
    size_t i = gp.getIndex(d);
    return ((i + 1) / 2) % 2 == 1;
  }

  void ZeroCrossingRefinementFunctor::getChild(const base::HashGridPoint& gp,
                                               size_t d, bool left,
                                               base::HashGridPoint& child) const {
    unsigned int i = gp.getIndex(d);
    unsigned int l = gp.getLevel(d);
    unsigned int new_l = l + 1;
    unsigned int new_i = i * 2;
    if (left) {
      new_i -= 1;
    } else {
      new_i += 1;
    }
    child.set(d, new_l, new_i);
  }

  void ZeroCrossingRefinementFunctor::getParent(const base::HashGridPoint& gp,
                                                size_t d,
                                                base::HashGridPoint& par) const {
    unsigned int i = gp.getIndex(d);
    unsigned int l = gp.getLevel(d);
    unsigned int new_l = l - 1;
    unsigned int new_i = (i + 1) / 2;
    if (new_i % 2 == 0) {
      new_i -= 1;
    }
    par.set(d, new_l, new_i);
  }

}  // namespace datadriven
}  // namespace sgpp
