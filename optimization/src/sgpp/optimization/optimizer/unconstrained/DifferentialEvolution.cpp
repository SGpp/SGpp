// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <sgpp/optimization/optimizer/unconstrained/DifferentialEvolution.hpp>

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

namespace sgpp {
namespace optimization {
namespace optimizer {

DifferentialEvolution::DifferentialEvolution(const base::ScalarFunction& f, size_t maxFcnEvalCount,
                                             size_t populationSize, double crossoverProbability,
                                             double scalingFactor, size_t idleGenerationsCount,
                                             double avgImprovementThreshold,
                                             double maxDistanceThreshold)
    : UnconstrainedOptimizer(f, nullptr, nullptr, maxFcnEvalCount),
      populationSize((populationSize > 0) ? populationSize : 10 * f.getNumberOfParameters()),
      crossoverProbability(crossoverProbability),
      scalingFactor(scalingFactor),
      idleGenerationsCount(idleGenerationsCount),
      avgImprovementThreshold(avgImprovementThreshold),
      maxDistanceThreshold(maxDistanceThreshold) {}

DifferentialEvolution::DifferentialEvolution(const DifferentialEvolution& other)
    : UnconstrainedOptimizer(other),
      populationSize(other.populationSize),
      crossoverProbability(other.crossoverProbability),
      scalingFactor(other.scalingFactor),
      idleGenerationsCount(other.idleGenerationsCount),
      avgImprovementThreshold(other.avgImprovementThreshold),
      maxDistanceThreshold(other.maxDistanceThreshold) {}

DifferentialEvolution::~DifferentialEvolution() {}

void DifferentialEvolution::optimize() {
  base::Printer::getInstance().printStatusBegin("Optimizing (differential evolution)...");

  const size_t d = f->getNumberOfParameters();

  xOpt.resize(0);
  fOpt = std::numeric_limits<double>::quiet_NaN();
  xHist.resize(0, d);
  fHist.resize(0);

  // vector of individuals
  std::vector<base::DataVector> x1(populationSize, base::DataVector(d, 0.0));
  // another vector for the new population
  std::vector<base::DataVector> x2(x1);

  // pointers for swapping both vectors at the end of each iterations
  std::vector<base::DataVector>* xOld = &x1;
  std::vector<base::DataVector>* xNew = &x2;

  // function values at the points of the populations
  // (no need to swap those)
  base::DataVector fx(populationSize);

  // initial pseudorandom points
  for (size_t i = 0; i < populationSize; i++) {
    for (size_t t = 0; t < d; t++) {
      (*xOld)[i][t] = base::RandomNumberGenerator::getInstance().getUniformRN();
    }

    fx[i] = f->eval((*xOld)[i]);
  }

  // smallest function value in the population
  double fCurrentOpt = std::numeric_limits<double>::infinity();
  // index of the point with value fOpt
  size_t xOptIndex = 0;

  // iteration number of the last iteration with significant improvement
  size_t lastNonidleK = 0;
  // average of all function values
  double avg = 0.0;
  // average in the previous round
  double lastAvg = 0.0;
  // number of iterations
  const size_t maxK = std::max(static_cast<size_t>(2), N / populationSize) - 1;

  std::vector<std::vector<size_t>> a(maxK, std::vector<size_t>(populationSize, 0)), b = a, c = a,
                                                                                    j = a;
  std::vector<std::vector<base::DataVector>> prob(
      maxK, std::vector<base::DataVector>(populationSize, base::DataVector(d, 0)));

  // pregenerate all pseudorandom numbers because the
  // real algorithm is parallelized
  // (for comparability of results, and maybe the
  // RandomNumberGenerator isn't thread-safe)
  for (size_t k = 0; k < maxK; k++) {
    for (size_t i = 0; i < populationSize; i++) {
      do {
        a[k][i] = base::RandomNumberGenerator::getInstance().getUniformIndexRN(populationSize);
      } while (a[k][i] == i);

      do {
        b[k][i] = base::RandomNumberGenerator::getInstance().getUniformIndexRN(populationSize);
      } while ((b[k][i] == i) || (b[k][i] == a[k][i]));

      do {
        c[k][i] = base::RandomNumberGenerator::getInstance().getUniformIndexRN(populationSize);
      } while ((c[k][i] == i) || (c[k][i] == a[k][i]) || (c[k][i] == b[k][i]));

      j[k][i] = base::RandomNumberGenerator::getInstance().getUniformIndexRN(d);

      for (size_t t = 0; t < d; t++) {
        if (t != j[k][i]) {
          prob[k][i][t] = base::RandomNumberGenerator::getInstance().getUniformRN();
        }
      }
    }
  }

  // "real" algorithm loop
  for (size_t k = 0; k < maxK; k++) {
    // abbreviations
    const std::vector<size_t>&a_k = a[k], &b_k = b[k], &c_k = c[k];
    const std::vector<size_t>& j_k = j[k];
    const std::vector<base::DataVector>& prob_k = prob[k];

#pragma omp parallel shared(k, a_k, b_k, c_k, j_k, prob_k, xOld, fx, fCurrentOpt, xOptIndex, xNew)
    {  // NOLINT(whitespace/braces)
      base::DataVector y(d);
      base::ScalarFunction* curFPtr = f.get();
#ifdef _OPENMP
      std::unique_ptr<base::ScalarFunction> curF;

      if (omp_get_max_threads() > 1) {
        f->clone(curF);
        curFPtr = curF.get();
      }

#endif /* _OPENMP */

// for each point in the population
#pragma omp for schedule(dynamic)

      for (size_t i = 0; i < populationSize; i++) {
        const size_t &cur_a = a_k[i], &cur_b = b_k[i], &cur_c = c_k[i];
        const size_t& cur_j = j_k[i];
        const base::DataVector& prob_ki = prob_k[i];
        bool inDomain = true;

        // for each dimension
        for (size_t t = 0; t < d; t++) {
          const double& curProb = prob_ki[t];

          if ((t == cur_j) || (curProb < crossoverProbability)) {
            // mutate point in this dimension
            y[t] = (*xOld)[cur_a][t] + scalingFactor * ((*xOld)[cur_b][t] - (*xOld)[cur_c][t]);
          } else {
            // don't mutate point in this dimension
            y[t] = (*xOld)[i][t];
          }

          // mutated point is out of bounds ==> discard
          if ((y[t] < 0.0) || (y[t] > 1.0)) {
            inDomain = false;
            break;
          }
        }

        // evaluate mutated point (if not out of bounds)
        const double fy = (inDomain ? curFPtr->eval(y) : std::numeric_limits<double>::infinity());

        if (fy < fx[i]) {
// function_value is better ==> replace point with mutated one
#pragma omp critical
          {
            fx[i] = fy;

            if (fy < fCurrentOpt) {
              xOptIndex = i;
              fCurrentOpt = fy;
            }
          }

          for (size_t t = 0; t < d; t++) {
            (*xNew)[i][t] = y[t];
          }
        } else {
          // function value not better ==> keep old point
          for (size_t t = 0; t < d; t++) {
            (*xNew)[i][t] = (*xOld)[i][t];
          }
        }
      }
    }

    // swap populations
    std::swap(xOld, xNew);
    avg = 0.0;

    // calculate average function value
    for (size_t i = 0; i < populationSize; i++) {
      avg += fx[i];
    }

    avg /= static_cast<double>(populationSize);

    if (lastAvg - avg >= avgImprovementThreshold) {
      // significant improvement
      lastNonidleK = k;
    } else if (k - lastNonidleK >= idleGenerationsCount) {
      // last significant improvement too long ago
      // ==> calculate maximum squared distance of all points
      // to the best one
      double maxDistance2 = 0.0;

      for (size_t i = 0; i < populationSize; i++) {
        if (i == xOptIndex) {
          continue;
        }

        double distance2 = 0.0;

        for (size_t t = 0; t < d; t++) {
          distance2 +=
              ((*xOld)[i][t] - (*xOld)[xOptIndex][t]) * ((*xOld)[i][t] - (*xOld)[xOptIndex][t]);
        }

        if (distance2 > maxDistance2) {
          maxDistance2 = distance2;
        }
      }

      // stopping criterion
      if (std::sqrt(maxDistance2) < maxDistanceThreshold) {
        break;
      }
    }

    // save average in last_avg
    lastAvg = avg;

    // status message
    if (k % 10 == 0) {
      base::Printer::getInstance().printStatusUpdate(
          std::to_string(k) + " steps, f(x) = " + std::to_string(fCurrentOpt));
    }

    xHist.appendRow((*xOld)[xOptIndex]);
    fHist.append(fCurrentOpt);
  }

  // optimal point
  xOpt.resize(d);
  xOpt = (*xOld)[xOptIndex];
  fOpt = fCurrentOpt;

  base::Printer::getInstance().printStatusUpdate(std::to_string(maxK) +
                                                 " steps, f(x) = " + std::to_string(fCurrentOpt));
  base::Printer::getInstance().printStatusEnd();
}

size_t DifferentialEvolution::getPopulationSize() const { return populationSize; }

void DifferentialEvolution::setPopulationSize(size_t populationSize) {
  this->populationSize = populationSize;
}

double DifferentialEvolution::getCrossoverProbability() const { return crossoverProbability; }

void DifferentialEvolution::setCrossoverProbability(double crossoverProbability) {
  this->crossoverProbability = crossoverProbability;
}

double DifferentialEvolution::getScalingFactor() const { return scalingFactor; }

void DifferentialEvolution::setScalingFactor(double scalingFactor) {
  this->scalingFactor = scalingFactor;
}

size_t DifferentialEvolution::getIdleGenerationsCount() const { return idleGenerationsCount; }

void DifferentialEvolution::setIdleGenerationsCount(size_t idleGenerationsCount) {
  this->idleGenerationsCount = idleGenerationsCount;
}

double DifferentialEvolution::getAvgImprovementThreshold() const {
  return avgImprovementThreshold;
}

void DifferentialEvolution::setAvgImprovementThreshold(double avgImprovementThreshold) {
  this->avgImprovementThreshold = avgImprovementThreshold;
}

double DifferentialEvolution::getMaxDistanceThreshold() const { return maxDistanceThreshold; }

void DifferentialEvolution::setMaxDistanceThreshold(double maxDistanceThreshold) {
  this->maxDistanceThreshold = maxDistanceThreshold;
}

void DifferentialEvolution::clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const {
  clone = std::unique_ptr<UnconstrainedOptimizer>(new DifferentialEvolution(*this));
}
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp
