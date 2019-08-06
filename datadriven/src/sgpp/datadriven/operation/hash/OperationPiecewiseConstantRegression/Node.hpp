// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <vector>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>

namespace sgpp {
namespace datadriven {
namespace PiecewiseConstantRegression {

class Node {
 private:
  std::vector<double> x;
  std::vector<double> h;

  size_t dim;

  std::vector<size_t> supportIndizes;

  std::unique_ptr<Node> leftChild;
  std::unique_ptr<Node> rightChild;
  size_t childDim;

  double surplus;

  base::DataMatrix& dataset;
  base::DataVector& values;

  size_t childCount;

  bool verbose;

 public:
  static uint64_t integratedNodes;

  static uint64_t hierarchizeMaxLevel;

  Node(std::vector<double> x, std::vector<double> h, std::vector<size_t> supportIndizes,
       base::DataMatrix& dataset, base::DataVector& values, bool verbose = false);

  std::vector<size_t> getSupportIndizes(std::vector<double>& x, std::vector<double>& h,
                                        std::vector<size_t>& parentSupport);

  double getAverage(std::vector<size_t>& supportIndizes, std::vector<double>& x,
                     std::vector<double>& h);

  double getMSE(std::vector<size_t>& supportIndizes, std::vector<double>& x,
                 std::vector<double>& h, double supportValue);

  std::unique_ptr<Node> hierarchizeChild(std::vector<double>& x, std::vector<double>& h,
                                         std::vector<size_t>& parentSupport, double supportValue,
                                         double targetMSE, size_t targetMaxLevel, size_t nextDim,
                                         size_t levelLimit);

  void hierarchize(double targetMSE, size_t targetMaxLevel, double parentValue = 0.0,
                   size_t refineDim = 0, size_t levelLimit = 0);

  std::vector<double> getLeftChildX(size_t d);

  std::vector<double> getRightChildX(size_t d);

  std::vector<double> getChildH(size_t d);

  double evaluate(std::vector<double>& point);

  double integrate(sgpp::base::GridPoint& gridPoint, size_t& integratedNodes,
                    size_t levelLimit = 0);

  uint64_t getChildCount();

  uint64_t getHierarchizationMaxLevel();

  double getMSE() {
    double mse = 0.0;

    for (size_t i = 0; i < dataset.getNrows(); i++) {
      std::vector<double> point;
      dataset.getRow(i, point);
      double eval = this->evaluate(point);
      mse += (eval - values[i]) * (eval - values[i]);
    }

    return mse;
  }
};
}  // namespace PiecewiseConstantRegression
}  // namespace datadriven
}  // namespace sgpp
