// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/tools/PolynomialChaosExpansion.hpp>

namespace sgpp {
namespace datadriven {
PolynomialChaosExpansion::PolynomialChaosExpansion(std::function<double(const base::DataVector&)> f,
                                                   int order,
                                                   sgpp::base::DistributionsVector distributions)
    : order(order),
      distributions(distributions),
      ranges(std::vector<std::pair<double, double>>(distributions.getSize())),
      alpha(base::DataVector(distributions.getSize(), 0)),
      beta(base::DataVector(distributions.getSize(), 0)) {
  types = std::vector<sgpp::datadriven::distributionType>(distributions.getSize());
  for (std::vector<distributionType>::size_type i = 0; i < types.size(); ++i) {
    auto characteristics = distributions.get(i)->getCharacteristics();
    if (distributions.get(i)->getType() == sgpp::base::DistributionType::Normal) {
      types[i] = sgpp::datadriven::distributionType::Normal;
      ranges[i].first = -9;
      ranges[i].second = 9;
      auto dist1 = std::make_shared<sgpp::base::DistributionNormal>(0, 1);
      standardvec.push_back(dist1);
    } else if (distributions.get(i)->getType() == sgpp::base::DistributionType::TruncNormal) {
      types[i] = sgpp::datadriven::distributionType::Normal;
      ranges[i].first = -9;
      ranges[i].second = 9;
      auto dist1 = std::make_shared<sgpp::base::DistributionNormal>(0, 1);
      standardvec.push_back(dist1);
    } else if (distributions.get(i)->getType() == sgpp::base::DistributionType::Lognormal) {
      types[i] = sgpp::datadriven::distributionType::Normal;
      ranges[i].first = -9;
      ranges[i].second = 9;
      auto dist1 = std::make_shared<sgpp::base::DistributionNormal>(0, 1);
      standardvec.push_back(dist1);
    } else if (distributions.get(i)->getType() == sgpp::base::DistributionType::Uniform) {
      types[i] = sgpp::datadriven::distributionType::Uniform;
      this->ranges[i].first = -(1.0);
      this->ranges[i].second = (1.0);
      auto dist1 = std::make_shared<sgpp::base::DistributionUniform>(-1, 1);
      standardvec.push_back(dist1);
    } else if (distributions.get(i)->getType() == sgpp::base::DistributionType::TruncGamma) {
      types[i] = sgpp::datadriven::distributionType::Gamma;
      this->alpha[i] = characteristics[0];
      this->ranges[i].first = 0;
      this->ranges[i].second = characteristics[1];
      auto dist1 = std::make_shared<sgpp::base::DistributionTruncGamma>(characteristics[0],
                                                                        characteristics[1]);
      standardvec.push_back(dist1);
    } else if (distributions.get(i)->getType() == sgpp::base::DistributionType::Beta) {
      types[i] = sgpp::datadriven::distributionType::Beta;
      this->ranges[i].first = -(1.0);
      this->ranges[i].second = (1.0);
      this->alpha[i] = characteristics[0];
      this->beta[i] = characteristics[1];
      auto dist1 =
          std::make_shared<sgpp::base::DistributionBeta>(characteristics[0], characteristics[1]);
      standardvec.push_back(dist1);
    } else if (distributions.get(i)->getType() == sgpp::base::DistributionType::TruncExponential) {
      types[i] = sgpp::datadriven::distributionType::Exponential;
      this->ranges[i].first = 0;
      this->ranges[i].second = characteristics[0];
      auto dist1 = std::make_shared<sgpp::base::DistributionTruncExponential>(characteristics[0],
                                                                              characteristics[1]);
      standardvec.push_back(dist1);
    }
  }

  /*
   * normalization
   */
  func = [f, this](const base::DataVector& vec) {
    base::DataVector temp(types.size());
    for (std::vector<distributionType>::size_type i = 0; i < types.size(); ++i) {
      auto characteristics = this->distributions.get(i)->getCharacteristics();
      if (this->distributions.get(i)->getType() == sgpp::base::DistributionType::Normal) {
        temp[i] = vec[i] * characteristics[1] + characteristics[0];
      } else if (this->distributions.get(i)->getType() ==
                 sgpp::base::DistributionType::TruncNormal) {
        temp[i] = vec[i] * characteristics[1] + characteristics[0];
      } else if (this->distributions.get(i)->getType() == sgpp::base::DistributionType::Lognormal) {
        temp[i] = exp((characteristics[0]) + vec[i] * (characteristics[1]));
      } else if (this->distributions.get(i)->getType() == sgpp::base::DistributionType::Uniform) {
        temp[i] = (vec[i] * (characteristics[1] - characteristics[0]) / 2) +
                  ((characteristics[1] + characteristics[0]) / 2);
      } else if (this->distributions.get(i)->getType() ==
                 sgpp::base::DistributionType::TruncGamma) {
        temp[i] = vec[i];
      } else if (this->distributions.get(i)->getType() == sgpp::base::DistributionType::Beta) {
        temp[i] = vec[i];
      } else if (this->distributions.get(i)->getType() ==
                 sgpp::base::DistributionType::TruncExponential) {
        temp[i] = vec[i];
      }
    }
    return f(temp);
  };
  // order:hermite,jacobi,legendre,laguerre,genlaguerre
  this->denoms = std::vector<std::function<double(double, size_t)>>{
      {[](double j, size_t i) { return std::tgamma(j + 1.0); }},
      {[this](double j, size_t i) {
        return (1.0 / (2.0 * j + this->alpha[i] + this->beta[i] + 1.0)) *
               ((std::tgamma(j + this->alpha[i] + 1.0) * std::tgamma(j + this->beta[i] + 1.0)) /
                (std::tgamma(j + this->alpha[i] + this->beta[i] + 1.0) * std::tgamma(j + 1.0))) /
               ((std::tgamma(this->alpha[i] + 1) * std::tgamma(this->beta[i] + 1)) /
                std::tgamma(this->alpha[i] + this->beta[i] + 2));
      }},
      {[](double j, size_t i) { return 1.0 / ((2.0 * j) + 1.0); }},
      {[](double j, size_t i) { return 1.0; }},
      {[this](double j, size_t i) {
        return (std::tgamma(j + this->alpha[i] + 1.0) /
                (std::tgamma(j + 1.0) * std::tgamma(this->alpha[i] + 1.0)));
      }}};
  this->evals = std::vector<std::function<double(double, double, size_t)>>{
      {[this](int n, double x, size_t i) { return evalHermite(n, x); }},
      {[this](int n, double x, size_t i) { return evalJacobi(n, x, i); }},
      {[this](int n, double x, size_t i) { return evalLegendre(n, x); }},
      {[this](int n, double x, size_t i) { return evalLaguerre(n, x); }},
      {[this](int n, double x, size_t i) { return evalGenLaguerre(n, x, i); }}};
}

PolynomialChaosExpansion::~PolynomialChaosExpansion() {}

double PolynomialChaosExpansion::evalLegendre(int n, double x) {
  if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return x;
  } else {
    double next = 0.0;
    double last = 1.0;
    double curr = x;
    for (double i = 1.0; i < n; ++i) {
      next = ((2.0 * i + 1.0) / (i + 1.0)) * x * curr - (i / (i + 1.0)) * last;
      last = curr;
      curr = next;
    }
    return next;
  }
}

double PolynomialChaosExpansion::evalHermite(int n, double x) {
  if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return x;
  } else {
    double next = 0.0;
    double last = 1.0;
    double curr = x;
    for (int i = 1.0; i < n; ++i) {
      next = x * curr - i * last;
      last = curr;
      curr = next;
    }
    return next;
  }
}

double PolynomialChaosExpansion::evalLaguerre(int n, double x) {
  if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return 1.0 - x;
  } else {
    double next = 0.0;
    double last = 1.0;
    double curr = 1.0 - x;
    for (double i = 1.0; i < n; ++i) {
      next = ((2.0 * i + 1.0 - x) * curr - i * last) / (i + 1.0);
      last = curr;
      curr = next;
    }
    return next;
  }
}

double PolynomialChaosExpansion::evalJacobi(int n, double x, size_t i) {
  if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return (alpha[i] + 1.0) + (alpha[i] + beta[i] + 2.0) * ((x - 1.0) / 2.0);
  } else {
    double next = 0.0;
    double last = 1.0;
    double curr = (alpha[i] + 1.0) + (alpha[i] + beta[i] + 2.0) * ((x - 1.0) / 2.0);
    for (size_t i = 2; i <= static_cast<size_t>(n); ++i) {
      double q1 = ((2.0 * static_cast<double>(i) + alpha[i] + beta[i] - 1.0) *
                   ((2.0 * static_cast<double>(i) + alpha[i] + beta[i]) * (2.0 * static_cast<double>(i) + alpha[i] + beta[i] - 2.0) * x +
                    std::pow(alpha[i], 2) - std::pow(beta[i], 2))) /
                  (2.0 * static_cast<double>(i) * (static_cast<double>(i) + alpha[i] + beta[i]) * (2.0 * static_cast<double>(i) + alpha[i] + beta[i] - 2.0));
      double q2 =
          (2.0 * (static_cast<double>(i) + alpha[i] - 1.0) * (static_cast<double>(i) + beta[i] - 1.0) * (2.0 * static_cast<double>(i) + alpha[i] + beta[i])) /
          (2.0 * static_cast<double>(i) * (static_cast<double>(i) + alpha[i] + beta[i]) * (2.0 * static_cast<double>(i) + alpha[i] + beta[i] - 2.0));
      next = q1 * curr - q2 * last;
      last = curr;
      curr = next;
    }
    return next;
  }
}

double PolynomialChaosExpansion::evalGenLaguerre(int n, double x, size_t i) {
  if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return 1.0 + alpha[i] - x;
  } else {
    double next = 0.0;
    double last = 1.0;
    double curr = 1.0 + alpha[i] - x;
    for (size_t i = 1; i < static_cast<size_t>(n); ++i) {
      next = ((2.0 * static_cast<double>(i) + 1.0 + alpha[i] - x) * curr - (static_cast<double>(i) + alpha[i]) * last) / (static_cast<double>(i) + 1.0);
      last = curr;
      curr = next;
    }
    return next;
  }
}

double PolynomialChaosExpansion::sparseGridQuadrature(
    const std::function<double(const base::DataVector&)>& funct, int dim, int n, size_t quadOrder) {
  auto numfunc = [&funct](const base::DataVector& input,
                          const std::vector<std::pair<double, double>>& ranges) {
    base::DataVector temp(input.size());
    for (base::DataVector::size_type i = 0; i < input.size(); ++i) {
      temp[i] = input[i] * (ranges[i].second - ranges[i].first) + ranges[i].first;
    }
    return funct(temp);
  };
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createNakBsplineExtendedGrid(dim, 3));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  int i = 0;
  while (gridStorage.getSize() < static_cast<size_t>(n)) {
    gridStorage.clear();
    grid->getGenerator().regular(i);
    ++i;
  }
  sgpp::base::DataVector evals(gridStorage.getSize());
  base::DataVector p(dim);
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    for (int j = 0; j < dim; ++j) {
      p[j] = gp.getStandardCoordinate(j);
    }
    evals[i] = numfunc(p, ranges);
  }
  bool succHierarch = false;

  try {
    std::unique_ptr<base::OperationHierarchisation>(
        sgpp::op_factory::createOperationHierarchisation(*grid))
        ->doHierarchisation(evals);
    succHierarch = true;
  } catch (...) {
    succHierarch = false;
  }

  base::DataVector coeffs(evals.getSize());
  if (!succHierarch) {
    sgpp::base::HierarchisationSLE hierSLE(*grid);
    sgpp::base::sle_solver::Eigen sleSolver;

    // solve linear system
    if (!sleSolver.solve(hierSLE, evals, coeffs)) {
      // return 1;
    }
    // overwrite evals
    evals = coeffs;
  }

  std::unique_ptr<sgpp::base::OperationWeightedQuadrature> opWQ(
      sgpp::op_factory::createOperationWeightedQuadrature(*grid, quadOrder));
  double res = opWQ->doWeightedQuadrature(coeffs, standardvec);
  return res;
}
double PolynomialChaosExpansion::adaptiveQuadratureWeighted(
    const std::function<double(const base::DataVector&)>& funct, int dim, size_t n,
    size_t quadOrder) {
  auto numfunc = [&funct](const base::DataVector& input,
                          const std::vector<std::pair<double, double>>& ranges) {
    base::DataVector temp(input.size());
    for (base::DataVector::size_type i = 0; i < ranges.size(); ++i) {
      temp[i] = (input[i]) * (ranges[i].second - ranges[i].first) + ranges[i].first;
    }
    return funct(temp);
  };

  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createNakBsplineExtendedGrid(dim, 3));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(2);

  base::DataVector coeffs(gridStorage.getSize());
  base::DataVector funEvals(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    base::GridPoint& gp = gridStorage.getPoint(i);
    base::DataVector vec(dim);
    for (int j = 0; j < dim; ++j) {
      vec[j] = gp.getStandardCoordinate(j);
    }
    funEvals[i] = numfunc(vec, ranges);
  }
  std::vector<size_t> addedPoints;
  /**
   * Refine adaptively until number of points is reached.
   */
  while (gridStorage.getSize() < n) {
    base::SurplusRefinementFunctor functor(coeffs, 10);
    grid->getGenerator().refine(functor, &addedPoints);

    coeffs.resize(gridStorage.getSize());
    funEvals.resize(gridStorage.getSize());

    for (size_t i = 0; i < addedPoints.size(); i++) {
      size_t seq = addedPoints[i];
      base::GridPoint& gp = gridStorage.getPoint(seq);
      base::DataVector vec(dim);
      for (int j = 0; j < dim; ++j) {
        vec[j] = gp.getStandardCoordinate(j);
      }
      funEvals[seq] = numfunc(vec, ranges);
    }

    coeffs.copyFrom(funEvals);

    // try hierarchisation
    bool succHierarch = false;

    try {
      std::unique_ptr<base::OperationHierarchisation>(
          sgpp::op_factory::createOperationHierarchisation(*grid))
          ->doHierarchisation(coeffs);
      succHierarch = true;
    } catch (...) {
      succHierarch = false;
    }

    if (!succHierarch) {
      sgpp::base::HierarchisationSLE hierSLE(*grid);
      sgpp::base::sle_solver::Eigen sleSolver;

      // solve linear system
      if (!sleSolver.solve(hierSLE, funEvals, coeffs)) {
        // return 1;
      }
    }
    addedPoints.clear();
  }
  std::unique_ptr<sgpp::base::OperationWeightedQuadrature> opWQ(
      sgpp::op_factory::createOperationWeightedQuadrature(*grid, quadOrder));
  double res = opWQ->doWeightedQuadrature(coeffs, standardvec);
  return res;
}

std::vector<std::vector<int>> PolynomialChaosExpansion::multiIndex(int dimension, int order) {
  std::vector<std::vector<int>> index(static_cast<int>(std::pow(order + 1, dimension)));
  std::vector<int> curr(dimension);
  for (size_t j = 0; j < index.size(); ++j) {
    index[j] = curr;
    curr[0]++;
    for (int i = 0; i < dimension; ++i) {
      if (curr[dimension - 1] > order) {
        break;
      }
      if (curr[i] > order) {
        curr[i] -= order + 1;
        curr[i + 1]++;
      }
    }
    if (curr[dimension - 1] > order) {
      break;
    }
  }
  index.erase(std::remove_if(index.begin(), index.end(),
                             [order](std::vector<int> ee) {
                               return std::accumulate(ee.begin(), ee.end(), 0) > order;
                             }),
              index.end());

  return index;
}
base::DataVector PolynomialChaosExpansion::calculateCoefficients(int n, bool use_adaptive) {
  auto index = multiIndex(static_cast<int>(types.size()), order);
  base::DataVector result(index.size());
  // calculate aj for each entry in the multiIndex
  for (std::vector<std::vector<int>>::size_type j = 0; j < index.size(); ++j) {
    auto numfunc = [this, &index, &j](const base::DataVector& vec) {
      double prd = 1;
      for (base::DataVector::size_type i = 0; i < vec.getSize(); ++i) {
        prd *= evals[static_cast<int>(types[i])](index[j][i], vec[i], i);
      }
      return prd;
    };
    auto intfunc = [this, &numfunc](const base::DataVector& vec) {
      return numfunc(vec) * func(vec);
    };
    double num = 1;
    if (!use_adaptive) {
      num = sparseGridQuadrature(intfunc, static_cast<int>(types.size()), n, 100);
    } else if (use_adaptive) {
      num = adaptiveQuadratureWeighted(intfunc, static_cast<int>(types.size()), n, 100);
    }
    // calculate denominator
    double denom = 1.0;
    for (std::vector<distributionType>::size_type i = 0; i < types.size(); ++i) {
      denom *= denoms[static_cast<int>(types[i])](index[j][i], i);
    }
    double aj = num / denom;
    result[j] = aj;
  }
  this->coefficients = result;
  return result;
}

base::DataVector PolynomialChaosExpansion::getCoefficients() { return coefficients; }

void PolynomialChaosExpansion::clearCoefficients() { this->coefficients.clear(); }

double PolynomialChaosExpansion::evalExpansion(const base::DataVector& xi, int n,
                                               bool use_adaptive) {
  if (coefficients.empty()) {
    calculateCoefficients(n, use_adaptive);
  }
  base::DataVector temp(xi);
  for (std::vector<distributionType>::size_type i = 0; i < types.size(); ++i) {
    auto characteristics = this->distributions.get(i)->getCharacteristics();
    if (this->distributions.get(i)->getType() == sgpp::base::DistributionType::Normal) {
      temp[i] = (xi[i] - characteristics[0]) / characteristics[1];
    } else if (this->distributions.get(i)->getType() == sgpp::base::DistributionType::TruncNormal) {
      temp[i] = (xi[i] - characteristics[0]) / characteristics[1];
    } else if (this->distributions.get(i)->getType() == sgpp::base::DistributionType::Lognormal) {
      temp[i] = (log(xi[i]) - (characteristics[0])) / (characteristics[1]);
    } else if (this->distributions.get(i)->getType() == sgpp::base::DistributionType::Uniform) {
      temp[i] = (xi[i] - ((characteristics[1] + characteristics[0]) / 2)) /
                ((characteristics[1] - characteristics[0]) / 2);
    } else if (this->distributions.get(i)->getType() == sgpp::base::DistributionType::TruncGamma) {
      temp[i] = xi[i];
    } else if (this->distributions.get(i)->getType() == sgpp::base::DistributionType::Beta) {
      temp[i] = xi[i];
    } else if (this->distributions.get(i)->getType() ==
               sgpp::base::DistributionType::TruncExponential) {
      temp[i] = xi[i];
    }
  }
  auto index = multiIndex(static_cast<int>(types.size()), order);
  double sum = 0.0;
  for (std::vector<std::vector<int>>::size_type j = 0; j < index.size(); ++j) {
    double prod = 1.0;
    for (std::vector<int>::size_type i = 0; i < index[j].size(); ++i) {
      prod *= evals[static_cast<int>(types[i])](index[j][i], temp[i], i);
    }
    sum += prod * coefficients[j];
  }
  return sum;
}

double PolynomialChaosExpansion::getL2Error(int n, bool use_adaptive) {
  auto gen = [this, &n, &use_adaptive]() {
    base::DataVector randvec = distributions.sample();
    base::DataVector transvec(randvec.getSize());
    for (std::vector<distributionType>::size_type i = 0; i < types.size(); ++i) {
      auto characteristics = this->distributions.get(i)->getCharacteristics();
      if (this->distributions.get(i)->getType() == sgpp::base::DistributionType::Normal) {
        transvec[i] = (randvec[i] - characteristics[0]) / characteristics[1];
      } else if (this->distributions.get(i)->getType() ==
                 sgpp::base::DistributionType::TruncNormal) {
        transvec[i] = (randvec[i] - characteristics[0]) / characteristics[1];
      } else if (this->distributions.get(i)->getType() == sgpp::base::DistributionType::Lognormal) {
        transvec[i] = (log(randvec[i]) - (characteristics[0])) / (characteristics[1]);
      } else if (this->distributions.get(i)->getType() == sgpp::base::DistributionType::Uniform) {
        transvec[i] = (randvec[i] - ((characteristics[1] + characteristics[0]) / 2)) /
                      ((characteristics[1] - characteristics[0]) / 2);
      } else if (this->distributions.get(i)->getType() ==
                 sgpp::base::DistributionType::TruncGamma) {
        transvec[i] = randvec[i];
      } else if (this->distributions.get(i)->getType() == sgpp::base::DistributionType::Beta) {
        transvec[i] = randvec[i];
      } else if (this->distributions.get(i)->getType() ==
                 sgpp::base::DistributionType::TruncExponential) {
        transvec[i] = randvec[i];
      }
    }
    return std::pow(func(transvec) - (evalExpansion(randvec, n, use_adaptive)), 2);
  };
  size_t num = 100000;
  base::DataVector results(num);
  std::generate(results.begin(), results.end(), gen);
  return results.sum() / static_cast<double>(num);
}
double PolynomialChaosExpansion::getMean(int n, bool use_adaptive) {
  if (coefficients.empty()) {
    calculateCoefficients(n, use_adaptive);
  }
  return coefficients[0];
}
double PolynomialChaosExpansion::getVariance(int n, bool use_adaptive) {
  if (coefficients.empty()) {
    calculateCoefficients(n, use_adaptive);
  }
  // order:hermite,jacobi,legendre,laguerre,genlaguerre
  auto normalization = std::vector<std::function<double(double, size_t)>>{
      {[](double j, size_t i) { return std::sqrt(std::tgamma(j + 1.0)); }},
      {[this](double j, size_t i) {
        return std::sqrt(
            (1.0 / (2.0 * j + this->alpha[i] + this->beta[i] + 1.0)) *
            ((std::tgamma(j + this->alpha[i] + 1.0) * std::tgamma(j + this->beta[i] + 1.0)) /
             (std::tgamma(j + this->alpha[i] + this->beta[i] + 1.0) * std::tgamma(j + 1.0))) /
            (((std::tgamma(this->alpha[i] + 1) * std::tgamma(this->beta[i] + 1)) /
              std::tgamma(this->alpha[i] + this->beta[i] + 2))));
      }},
      {[](double j, size_t i) { return std::sqrt(1.0 / ((2.0 * j) + 1.0)); }},
      {[](double j, size_t i) { return 1.0; }},
      {[this](double j, size_t i) {
        return std::sqrt(std::tgamma(j + this->alpha[i] + 1.0) /
                         (std::tgamma(j + 1.0) * std::tgamma(this->alpha[i] + 1.0)));
      }}};
  base::DataVector temp(coefficients.getSize());
  temp.copyFrom(coefficients);
  temp.set(0, 0.0);
  auto index = multiIndex(static_cast<int>(types.size()), order);
  for (std::vector<std::vector<int>>::size_type j = 0; j < index.size(); ++j) {
    for (std::vector<int>::size_type i = 0; i < index[j].size(); ++i) {
      temp[j] *= normalization[static_cast<int>(types[i])](index[j][i], i);
    }
  }
  temp.sqr();
  return temp.sum();
}
}  // namespace datadriven
}  // namespace sgpp
