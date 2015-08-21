// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/unconstrained/DifferentialEvolution.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      DifferentialEvolution::DifferentialEvolution(ObjectiveFunction& f,
          size_t maxFcnEvalCount, size_t populationSize,
          float_t crossoverProbability, float_t scalingFactor,
          size_t idleGenerationsCount, float_t avgImprovementThreshold,
          float_t maxDistanceThreshold) :
        UnconstrainedOptimizer(f, maxFcnEvalCount),
        populationSize((populationSize > 0) ? populationSize :
                       10 * f.getDimension()),
        crossoverProbability(crossoverProbability),
        scalingFactor(scalingFactor),
        idleGenerationsCount(idleGenerationsCount),
        avgImprovementThreshold(avgImprovementThreshold),
        maxDistanceThreshold(maxDistanceThreshold) {
      }

      float_t DifferentialEvolution::optimize(base::DataVector& xOpt) {
        printer.printStatusBegin("Optimizing (differential evolution)...");

        const size_t d = f.getDimension();

        // vector of individuals
        std::vector<base::DataVector> x1(populationSize,
                                         base::DataVector(d, 0.0));
        // another vector for the new population
        std::vector<base::DataVector> x2(x1);

        // pointers for swapping both vectors at the end of each iterations
        std::vector<base::DataVector>* xOld = &x1;
        std::vector<base::DataVector>* xNew = &x2;

        // function values at the points of the populations
        // (no need to swape those)
        base::DataVector fx(populationSize);

        // initial pseudorandom points
        for (size_t i = 0; i < populationSize; i++) {
          for (size_t t = 0; t < d; t++) {
            (*xOld)[i][t] = randomNumberGenerator.getUniformRN();
          }

          fx[i] = f.eval((*xOld)[i]);
        }

        // smallest function value in the population
        float_t fOpt = INFINITY;
        // index of the point with value fOpt
        size_t xOptIndex = 0;
        // iteration number of the last iteration with significant improvement
        size_t lastNonidleK = 0;
        // average of all function values
        float_t avg = 0.0;
        // average in the previous round
        float_t lastAvg = 0.0;
        // number of iterations
        size_t maxK = std::max(static_cast<size_t>(2), N / populationSize) - 1;

        std::vector<std::vector<size_t>> a(maxK, std::vector<size_t>(
                                             populationSize, 0)),
                                                             b = a, c = a, j = a;
        std::vector<std::vector<base::DataVector>> prob(maxK,
                                                std::vector<base::DataVector>(populationSize,
                                                    base::DataVector(d, 0)));

        // pregenerate all pseudorandom numbers because the
        // real algorithm is parallelized
        // (for comparability of results, and maybe the
        // RandomNumberGenerator isn't thread-safe)
        for (size_t k = 0; k < maxK; k++) {
          for (size_t i = 0; i < populationSize; i++) {
            do {
              a[k][i] = randomNumberGenerator.
                        getUniformIndexRN(populationSize);
            } while (a[k][i] == i);

            do {
              b[k][i] = randomNumberGenerator.
                        getUniformIndexRN(populationSize);
            } while ((b[k][i] == i) || (b[k][i] == a[k][i]));

            do {
              c[k][i] = randomNumberGenerator.
                        getUniformIndexRN(populationSize);
            } while ((c[k][i] == i) || (c[k][i] == a[k][i]) ||
                     (c[k][i] == b[k][i]));

            j[k][i] = randomNumberGenerator.getUniformIndexRN(d);

            for (size_t t = 0; t < d; t++) {
              if (t != j[k][i]) {
                prob[k][i][t] = randomNumberGenerator.getUniformRN();
              }
            }
          }
        }

        // "real" algorithm loop
        for (size_t k = 0; k < maxK; k++) {
          // abbreviations
          const std::vector<size_t>& a_k = a[k], &b_k = b[k], &c_k = c[k];
          const std::vector<size_t>& j_k = j[k];
          const std::vector<base::DataVector>& prob_k = prob[k];

          #pragma omp parallel shared(k, a_k, b_k, c_k, j_k, prob_k, \
          xOld, fx, fOpt, xOptIndex, xNew) default(none)
          {
            base::DataVector y(d);
            ObjectiveFunction* curFPtr = &f;
#ifdef _OPENMP
            std::unique_ptr<ObjectiveFunction> curF;

            if (omp_get_max_threads() > 1) {
              f.clone(curF);
              curFPtr = curF.get();
            }

#endif /* _OPENMP */

            // for each point in the population
            #pragma omp for schedule(dynamic)

            for (size_t i = 0; i < populationSize; i++) {
              const size_t& cur_a = a_k[i], &cur_b = b_k[i], &cur_c = c_k[i];
              const size_t& cur_j = j_k[i];
              const base::DataVector& prob_ki = prob_k[i];
              bool inDomain = true;

              // for each dimension
              for (size_t t = 0; t < d; t++) {
                const float_t& curProb = prob_ki[t];

                if ((t == cur_j) || (curProb < crossoverProbability)) {
                  // mutate point in this dimension
                  y[t] = (*xOld)[cur_a][t] +
                         scalingFactor * ((*xOld)[cur_b][t] -
                                          (*xOld)[cur_c][t]);
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
              const float_t fy = (inDomain ? curFPtr->eval(y) : INFINITY);

              if (fy < fx[i]) {
                // function_value is better ==> replace point with mutated one
                #pragma omp critical
                {
                  fx[i] = fy;

                  if (fy < fOpt) {
                    xOptIndex = i;
                    fOpt = fy;
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

          avg /= static_cast<float_t>(populationSize);

          if (lastAvg - avg >= avgImprovementThreshold) {
          // significant improvement
          lastNonidleK = k;
        } else if (k - lastNonidleK >= idleGenerationsCount) {
          // last significant improvement too long ago
          // ==> calculate maximum squared distance of all points
          // to the best one
          float_t maxDistance2 = 0.0;

          for (size_t i = 0; i < populationSize; i++) {
              if (i == xOptIndex) {
                continue;
              }

              float_t distance2 = 0.0;

              for (size_t t = 0; t < d; t++) {
                distance2 += ((*xOld)[i][t] - (*xOld)[xOptIndex][t]) *
                             ((*xOld)[i][t] - (*xOld)[xOptIndex][t]);
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
          printer.printStatusUpdate(std::to_string(k) + " steps, f(x) = " +
                                    std::to_string(fOpt));
          }
        }

        // optimal point
        xOpt.resize(d);
        xOpt = (*xOld)[xOptIndex];

        printer.printStatusUpdate(std::to_string(maxK) + " steps, f(x) = " +
                                  std::to_string(fOpt));
        printer.printStatusEnd();

        return fOpt;
      }

      size_t DifferentialEvolution::getPopulationSize() const {
        return populationSize;
      }

      void DifferentialEvolution::setPopulationSize(size_t populationSize) {
        this->populationSize = populationSize;
      }

    }
  }
}
