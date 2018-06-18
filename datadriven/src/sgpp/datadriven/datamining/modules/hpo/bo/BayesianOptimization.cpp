/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * BayesianOptimization.cpp
 *
 *  Created on: Feb 2, 2018
 *      Author: Eric Koepke
 */

#include <sgpp/optimization/sle/solver/Eigen.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/bo/BayesianOptimization.hpp>
#include <iostream>
#include <sgpp/optimization/sle/solver/BiCGStab.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/sle/solver/GaussianElimination.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/optimizer/unconstrained/MultiStart.hpp>

namespace sgpp {
namespace datadriven {

BayesianOptimization::BayesianOptimization(const std::vector<BOConfig> &initialConfigs)
    : kernelmatrix(initialConfigs.size(), initialConfigs.size()),
      gleft(),
      transformedOutput(),
      rawScores(initialConfigs.size()),
      screwedvar(false),
      maxofmax(0),
      allConfigs(initialConfigs) {
  scales = base::DataVector(initialConfigs.front().getNPar() + 1, 1);
  for (size_t i = 0; i < allConfigs.size(); ++i) {
    rawScores[i] = transformScore(allConfigs[i].getScore());
  }
  rawScores.normalize();
  rawScores.sub(base::DataVector(rawScores.size(),
                                 rawScores.sum() / static_cast<double>(rawScores.size())));
  bestsofar = rawScores.min();
  fitScales();
}

double BayesianOptimization::transformScore(double original) {
  return -1 / (1 + original);
}

double BayesianOptimization::kernel(double distance) {
  double kernelwidth = 8.0;   // offset for dimensional scaling
  return exp(-distance * kernelwidth / 2);
}

double BayesianOptimization::acquisitionEI(base::DataVector knew, double kself, double bestsofar) {
  double mean = knew.dotProduct(transformedOutput);

  base::DataVector tmp(knew);
  solveCholeskySystem(gleft, tmp);

  base::DataVector check(tmp.size());
  kernelmatrix.mult(tmp, check);
  check.sub(knew);
  double max = check.maxNorm();
  maxofmax = std::fmax(max, maxofmax);
  double var = kself - knew.dotProduct(tmp);
  if (var > 1 || var < 0) {
    // std::cout << "Error: wrong variance: "<<var;
    // <<", knew:"<<knew.toString()<<", temp: "<<tmp.toString()<<std::endl;
    screwedvar = true;
    return 0;
  }

  double z = (mean - (bestsofar - 0.001)) / var;
  return ((mean - (bestsofar - 0.001)) * (0.5 + 0.5 * std::erf(-z / 1.41))
      - var * 0.4 * std::exp(-0.5 * z * z));   // EDIT: is this calculated properly?
}

BOConfig *BayesianOptimization::main(BOConfig &prototype) {
  BOConfig nextconfig(prototype);
  optimization::WrapperScalarFunction wrapper(prototype.getContSize(),
                                              std::bind(&BayesianOptimization::acquisitionOuter,
                                                        this,
                                                        std::placeholders::_1));
  double min = std::numeric_limits<double>::infinity();
  BOConfig bestConfig;
  do {
    for (auto &config : allConfigs) {
      config.calcDiscDistance(nextconfig, scales);
    }
    optimization::optimizer::MultiStart optimizer(wrapper);
    optimizer.optimize();
    double optv = optimizer.getOptimalValue();
    if (prototype.getContSize() == 0) {
      optv = acquisitionOuter(base::DataVector());
    }
    if (optv < min) {
      min = optv;
      bestConfig = BOConfig(nextconfig);
      bestConfig.setCont(optimizer.getOptimalPoint());
    }
  } while (nextconfig.nextDisc());
  // std::cout << "Acquistion: " << min << std::endl;
  allConfigs.push_back(bestConfig);
  return &allConfigs.back();
}

double BayesianOptimization::acquisitionOuter(const base::DataVector &inp) {
  base::DataVector kernelrow(allConfigs.size());
  for (size_t i = 0; i < allConfigs.size(); i++) {
    kernelrow[i] = kernel(allConfigs[i].getTotalDistance(inp, scales));   // divided by 2
    // EDIT: removed kernel = 1 check, okay?
  }
  return acquisitionEI(kernelrow, 1, bestsofar);  // EDIT: kself + noise?
}

void BayesianOptimization::fitScales() {
  optimization::WrapperScalarFunction wrapper(scales.size(),
                                              std::bind(&BayesianOptimization::likelihood,
                                                        this,
                                                        std::placeholders::_1));
  // adjust resource allocation for optimizer here
  optimization::optimizer::MultiStart optimizer(wrapper, 2000, 200);
  optimizer.optimize();
  // std::cout << optimizer.getOptimalPoint().toString() << std::endl;
  // std::cout << optimizer.getOptimalValue() << std::endl;
  base::DataVector tmp(optimizer.getOptimalPoint());
  tmp.mult(0.1);   // factor for smooth updating of GP
  scales.mult(0.9);
  scales.add(tmp);
  // std::cout << scales.toString() << std::endl;
  double noise = pow(10, -scales.back() * 10);
  for (size_t i = 0; i < allConfigs.size(); ++i) {
    for (size_t k = 0; k < i; ++k) {
      double tmp = kernel(allConfigs[i].getScaledDistance(allConfigs[k], scales));
      kernelmatrix.set(k, i, tmp);
      kernelmatrix.set(i, k, tmp);
    }
    kernelmatrix.set(i, i, 1 + noise);
  }
  decomposeCholesky(kernelmatrix, gleft);
  transformedOutput = base::DataVector(rawScores);
  solveCholeskySystem(gleft, transformedOutput);
}

double BayesianOptimization::likelihood(const base::DataVector &inp) {
  double noise = pow(10, -inp.back() * 10);
  base::DataMatrix km(allConfigs.size(), allConfigs.size());
  for (size_t i = 0; i < allConfigs.size(); ++i) {
    for (size_t k = 0; k < i; ++k) {
      double tmp = kernel(allConfigs[i].getScaledDistance(allConfigs[k], inp));
      km.set(k, i, tmp);
      km.set(i, k, tmp);
    }
    km.set(i, i, 1 + noise);
  }
  base::DataMatrix gnew;
  decomposeCholesky(km, gnew);

  base::DataVector transformed(rawScores);
  solveCholeskySystem(gnew, transformed);
  double tmp = 0;
  for (size_t i = 0; i < allConfigs.size(); ++i) {
    tmp += std::log(gnew.get(i, i));
  }
  return 2 * tmp + rawScores.dotProduct(transformed);   // EDIT: correct likelihood?
}

void BayesianOptimization::updateGP() {
  double noise = pow(10, -scales.back() * 10);
  size_t size = kernelmatrix.getNcols();
  kernelmatrix.appendRow();
  kernelmatrix.appendCol(base::DataVector(size + 1));
  for (size_t i = 0; i < size; ++i) {
    double tmp = kernel(allConfigs[i].getScaledDistance(allConfigs.back(), scales));
    kernelmatrix.set(size, i, tmp);
    kernelmatrix.set(i, size, tmp);
    rawScores[i] = transformScore(allConfigs[i].getScore());
  }
  kernelmatrix.set(size, size, 1 + noise);
  rawScores.push_back(transformScore(allConfigs.back().getScore()));

  decomposeCholesky(kernelmatrix, gleft);

  rawScores.normalize();
  rawScores.sub(base::DataVector(rawScores.size(),
                                 rawScores.sum() / static_cast<double>(rawScores.size())));
  bestsofar = rawScores.min();
  transformedOutput = base::DataVector(rawScores);
  solveCholeskySystem(gleft, transformedOutput);


  // std::cout << "Maxofmax: " << maxofmax << std::endl;
  base::DataVector check2(transformedOutput.size());
  kernelmatrix.mult(transformedOutput, check2);
  check2.sub(rawScores);
  double max = check2.maxNorm();
  // std::cout << "Max: " << max << std::endl;
  // std::cout<<transformedOutput.toString()<<std::endl;

  // std::cout << "Var Screwed: " << screwedvar << std::endl;
  if (decomFailed || screwedvar || maxofmax > 0.1 || max > 0.1) {
    std::cout << "Numerical instabilities occured. This could lead to bad sampling.";
  }
  maxofmax = 0;
  screwedvar = false;
  decomFailed = false;
}

void BayesianOptimization::decomposeCholesky(base::DataMatrix &km, base::DataMatrix &gnew) {
  size_t n = km.getNrows();
  gnew = base::DataMatrix(n, n, 0);

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j <= i; j++) {
      double sum = km.get(i, j);
      for (size_t k = 0; k < j; k++) {
        sum = sum - gnew.get(i, k) * gnew.get(j, k);
      }
      if (i > j) {
        double value = sum / gnew.get(j, j);
        gnew.set(i, j, value);
      } else if (sum > 0) {
        double value = std::sqrt(sum);
        gnew.set(i, i, value);
      } else {
        decomFailed = true;
        gnew.set(i, i, 10e-8);
      }
    }
  }
}

void BayesianOptimization::solveCholeskySystem(base::DataMatrix &gmatrix, base::DataVector &x) {
  for (size_t i = 0; i < x.size(); i++) {
    x[i] = x[i] / gmatrix.get(i, i);
    for (size_t k = i + 1; k < x.size(); k++) {
      x[k] = x[k] - gmatrix.get(k, i) * x[i];
    }
  }

  for (int i = static_cast<int>(x.size()) - 1; i >= 0; i--) {
    x[i] = x[i] / gmatrix.get(static_cast<size_t>(i), static_cast<size_t>(i));
    for (int k = i - 1; k >= 0; k--) {
      x[k] = x[k] - gmatrix.get(static_cast<size_t>(i), static_cast<size_t>(k)) * x[i];
    }
  }
}
} /* namespace datadriven */
} /* namespace sgpp */
