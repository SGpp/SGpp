// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>


namespace sgpp {
namespace datadriven {

  DataBasedRefinementFunctor::
  DataBasedRefinementFunctor(std::vector<base::Grid*> grids,
                             std::vector<base::DataVector*> alphas,
                             std::vector<double> priors,
                             base::DataMatrix* data,
                             base::DataVector* targets,
                             size_t refinements_num,
                             bool levelPen,
                             std::vector<double> coeff_a,
                             double thresh) :
    grids(grids), alphas(alphas), priors(priors), evals(0, 0),
    data(data), targets(targets), h(grids.size()),
    means(), coeff_a(coeff_a), current_grid_index(0),
    refinements_num(refinements_num),
    threshold(thresh), level_penalize(levelPen) {
    // Default coeff_a = 1.0
    if (coeff_a.size() == 0) {
      this->coeff_a = std::vector<double>(grids.size(), 1.0);
    }

    // Compute H if data was provided
    if (data != nullptr && targets != nullptr) {
      computeH();
    }
  }

  double DataBasedRefinementFunctor::operator()(base::GridStorage& storage,
                                                size_t seq) const {
    size_t accum = 0;
    base::HashGridPoint& gp = storage.getPoint(seq);
    base::DataVector p(data->getNcols());

    // How many data points of H lie in the support of seq?
    for (size_t j = 0; j < h.at(current_grid_index).getNrows(); j++) {
      h.at(current_grid_index).getRow(j, p);
      if (isWithinSupport(gp, p)) {
        accum++;
      }
    }
    double score = static_cast<double>(accum);
    double levelSum = storage.getPoint(seq).getLevelSum();
    double levelW = pow(2.0, -levelSum);
    if (level_penalize) {
      score *= levelW;
    }
    return score;
  }

  double DataBasedRefinementFunctor::start() const {
    return 0.0;
  }

  size_t DataBasedRefinementFunctor::getRefinementsNum() const {
    return this->refinements_num;
  }

  double DataBasedRefinementFunctor::getRefinementThreshold() const {
    return this->threshold;
  }

  void DataBasedRefinementFunctor::setGridIndex(size_t grid_index) {
    this->current_grid_index = grid_index;
  }

  size_t DataBasedRefinementFunctor::getNumGrids() {
    return this->grids.size();
  }

  void DataBasedRefinementFunctor::setData(base::DataMatrix* data,
                                           base::DataVector* targets) {
    this->data = data;
    this->targets = targets;
  }

  void DataBasedRefinementFunctor::computeH() {
    // Evaluate all grids at all data points
    base::DataVector evalVec(data->getNrows());
    evals.resize(data->getNrows(), grids.size());
    for (size_t i = 0; i < grids.size(); i++) {
      std::unique_ptr<base::OperationMultipleEval>
        opEval(op_factory::createOperationMultipleEval(*grids.at(i),
                                                       *data));
      opEval->eval(*alphas.at(i), evalVec);
      evals.setColumn(i, evalVec);
      means.push_back(evalVec.sum() *
                      (priors.at(i) / static_cast<double>(data->getNrows())));
    }

    // Compute the sets H_k by pairwise H_kl for all class combiniations
    // of k != l
    for (size_t i = 0; i < grids.size(); i++) {
      h.push_back(base::DataMatrix());
      h.at(i).resize(0, data->getNcols());
      for (size_t j = 0; j < grids.size(); j++) {
        if (i == j) {
          continue;
        }
        computeHkl(h.at(i), i, j);
      }
    }
  }

  void DataBasedRefinementFunctor::computeHkl(base::DataMatrix& inters,
                                              size_t cl_ind1,
                                              size_t cl_ind2) {
    double val_1 = 0.0; double val_2 = 0.0;
    base::DataVector p(data->getNcols());
    for (size_t i = 0; i < evals.getNrows(); i++) {
      val_1 = evals.get(i, cl_ind1);
      val_2 = evals.get(i, cl_ind2);
      // If both PDFs surpass the threshold: mu * coeff_a, added the data
      // point
      if (val_1 > means.at(cl_ind1) * coeff_a.at(cl_ind1) &&
         val_2 > means.at(cl_ind2) * coeff_a.at(cl_ind2)) {
        data->getRow(i, p);
        inters.appendRow(p);
      }
    }
  }

  bool DataBasedRefinementFunctor::isWithinSupport(base::HashGridPoint& gp,
                                                   base::DataVector& point)
    const {
    for (size_t d = 0; d < point.getSize(); d++) {
      double coord = gp.getStandardCoordinate(d);
      size_t level = gp.getLevel(d);
      double step = 1.0 / pow(2.0, static_cast<double>(level));
      double lower = coord - step;
      double upper = coord + step;
      if (point.get(d) < lower || point.get(d) > upper) {
        return false;
      }
    }
    return true;
  }

  base::DataMatrix& DataBasedRefinementFunctor::getHk(size_t index) {
    return h.at(index);
  }

}  // namespace datadriven
}  // namespace sgpp
