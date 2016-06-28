// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/functors/classification/GridPointBasedRefinementFunctor.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>


namespace sgpp {
namespace datadriven {

  GridPointBasedRefinementFunctor::
  GridPointBasedRefinementFunctor(std::vector<base::Grid*> grids,
                                  std::vector<base::DataVector*> alphas,
                                  size_t r_num,
                                  bool level_penalize,
                                  bool pre_compute,
                                  double thresh) :
    grids(grids), alphas(alphas), current_grid_index(-1),
    refinements_num(r_num), threshold(thresh),
    level_penalize(level_penalize),
    pre_compute(pre_compute),
    pre_comp_evals() {
    for (size_t i = 0; i < grids.size(); i++) {
      pre_comp_evals.push_back(std::map<std::string, double>());
    }
  }

  double
  GridPointBasedRefinementFunctor::operator()(base::GridStorage& storage,
                                              size_t seq) const {
    double levelSum = storage.getPoint(seq).getLevelSum();
    double levelW = pow(2, -levelSum);
    double maxScore = 0;

    // Evaluations of all grids at (param) seq
    std::vector<double> gridEvals;

    // Get the evaluations of seq at all GridSave
    base::DataVector p(storage.getDimension());
    storage.getPoint(seq).getStandardCoordinates(p);
    if (pre_compute) {
      for (size_t i = 0; i < grids.size(); i++) {
        std::string key = p.toString();
        gridEvals.push_back(pre_comp_evals.at(i).at(key));
      }
    } else {
      for (size_t i = 0; i < grids.size(); i++) {
        std::unique_ptr<base::OperationEval>
          opEval(op_factory::createOperationEval(*grids.at(i)));
        gridEvals.push_back(opEval->eval(*alphas.at(i), p));
      }
    }

    // Do the scoring
    for (size_t i = 0; i < grids.size(); i++) {
      if (i == current_grid_index) {
        continue;
      }
      double diff = fabs(gridEvals.at(i) -
                         gridEvals.at(current_grid_index));
      double sqSum = pow(gridEvals.at(i), 2.0) +
        pow(gridEvals.at(current_grid_index), 2.0);
      double score = (sqSum / (1 + diff * diff));
      if (level_penalize) {
        score *= levelW;
      }
      if (score > maxScore) {
        maxScore = score;
      }
    }
    return maxScore;
  }

  void GridPointBasedRefinementFunctor::preComputeEvaluations() {
    base::DataVector p(grids.at(0)->getDimension());
    std::string key = "";
    double v = 0;
    for (size_t i = 0; i < grids.size(); i++) {
      std::unique_ptr<base::OperationEval>
        opEval(op_factory::createOperationEval(*grids.at(i)));
      pre_comp_evals.at(i).clear();

      for (size_t j = 0; j < grids.size(); j++) {
        for (size_t k = 0; k < grids.at(j)->getSize(); k++) {
          grids.at(j)->getStorage().getPoint(k).getStandardCoordinates(p);
          key = p.toString();
          if (pre_comp_evals.at(i).count(key) == 0) {
            v = opEval->eval(*alphas.at(i), p);
            pre_comp_evals.at(i).insert(std::pair<std::string, double>(key,
                                                                       v));
          }
        }
      }
    }
  }

  double GridPointBasedRefinementFunctor::start() const {
    return 0.0;
  }

  size_t GridPointBasedRefinementFunctor::getRefinementsNum() const {
    return this->refinements_num;
  }

  double GridPointBasedRefinementFunctor::getRefinementThreshold() const {
    return this->threshold;
  }

  void GridPointBasedRefinementFunctor::setGridIndex(size_t grid_index) {
    this->current_grid_index = grid_index;
  }

  size_t GridPointBasedRefinementFunctor::getNumGrids() {
    return this->grids.size();
  }

}
}
