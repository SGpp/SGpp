// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GROUPLASSOFUNCTION_HPP
#define GROUPLASSOFUNCTION_HPP

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/solver/sle/fista/RegularizationFunction.hpp>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sgpp {
namespace solver {

class GroupLassoFunction : public RegularizationFunction {
 public:
  GroupLassoFunction(double lambda, sgpp::base::GridStorage* gridStorage)
      : lambda(lambda), gridStorage(gridStorage) {}

  double eval(sgpp::base::DataVector weights) override {
    if (lastWeightSize != weights.getSize()) {
      calculateGroupIndices(weights);
      lastWeightSize = weights.getSize();
    }
    std::vector<double> norms;
    std::vector<size_t> groupSizes;
    std::tie(norms, groupSizes) = calculateNorms(weights);

    double penalty = 0.0;
    for (size_t i = 0; i < norms.size(); ++i) {
      penalty += std::sqrt(groupSizes[i]) * norms[i];
    }
    return penalty;
  }

  base::DataVector prox(const sgpp::base::DataVector& weights, double stepsize) override {
    if (lastWeightSize != weights.getSize()) {
      calculateGroupIndices(weights);
      lastWeightSize = weights.getSize();
    }
    std::vector<double> norms;
    std::vector<size_t> groupSizes;
    std::tie(norms, groupSizes) = calculateNorms(weights);

    base::DataVector proxVec = sgpp::base::DataVector(weights.getSize());
    for (size_t i = 0; i < proxVec.getSize(); ++i) {
      const auto curGroup = groupsVec[i];
      const auto curNorm = norms[curGroup];
      const auto curSize = groupSizes[curGroup];

      // Here we normalize with the root of the groupSize.
      // If we would not do that, we would drop larger groups with a higher probability than smaller
      // ones!
      // TODO(krenzls): Check formula!
      const double ss = lambda * stepsize * std::sqrt(curSize);
      const double multiplicator = std::max(1 - (ss) / curNorm, 0.0);
      proxVec[i] = multiplicator * weights[i];
    }
    return proxVec;
  }

 private:
  const double lambda;
  sgpp::base::GridStorage* gridStorage;
  std::unordered_map<std::vector<bool>, size_t> groups;
  std::vector<size_t> groupsVec;
  size_t lastWeightSize = 0;

  std::pair<std::vector<double>, std::vector<size_t>> calculateNorms(
      const sgpp::base::DataVector& weights) {
    auto groupSizes = std::vector<size_t>(groups.size());
    auto norms = std::vector<double>(groups.size());
    // First sum up the values.
    for (size_t i = 0; i < weights.getSize(); ++i) {
      const auto curGroup = groupsVec[i];
      groupSizes[curGroup]++;
      norms[curGroup] += weights[i] * weights[i];
    }
    // Then normalize.
    for (size_t i = 0; i < norms.size(); ++i) {
      norms[i] = std::sqrt(norms[i]);
    }

    return std::make_pair<std::vector<double>, std::vector<size_t>>(std::move(norms),
                                                                    std::move(groupSizes));
  }

  void calculateGroupIndices(const sgpp::base::DataVector& weights) {
    const size_t dim = gridStorage->getDimension();
    // Mapping between interaction-type and group number.
    // Note that there is a specialisation for std::vector<bool> that is implemented as a bitfield.
    groups = std::unordered_map<std::vector<bool>, size_t>();
    // Stores the group number of each weight.
    groupsVec = std::vector<size_t>(weights.getSize());
    size_t curGroupNum = 0;
    for (size_t i = 0; i < gridStorage->getSize(); ++i) {
      const auto p = gridStorage->get(i);
      auto coords = sgpp::base::DataVector(dim);
      p->getCoords(coords);

      // Check which coordinates are used, i.e. not equal to 0.5.
      auto dimsUsed = std::vector<bool>(dim);
      for (size_t j = 0; j < dim; ++j) {
        dimsUsed[j] = coords[j] != 0.5;
      }
      auto it = groups.find(dimsUsed);
      if (it == groups.end()) {
        // New group discovered, find a number for it.
        groups.emplace(dimsUsed, curGroupNum);
        groupsVec[i] = curGroupNum;
        ++curGroupNum;
      } else {
        groupsVec[i] = groups[dimsUsed];
      }
    }
  }
};

}  //  namespace solver
}  //  namepsace sgpp

#endif  // GROUPLASSOFUNCTION_HPP
