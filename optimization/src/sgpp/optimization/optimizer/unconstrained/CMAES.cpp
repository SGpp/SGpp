// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/unconstrained/CMAES.hpp>
#include <sgpp/optimization/tools/Math.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <algorithm>
#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      CMAES::CMAES(ObjectiveFunction& f,
                   size_t maxFcnEvalCount,
                   size_t populationSize,
                   size_t recombinationSize,
                   float_t initialStepSize) :
        UnconstrainedOptimizer(f, maxFcnEvalCount),
        lambda((populationSize > 0) ? populationSize :
               4 + static_cast<size_t>(3.0 * std::log(f.getDimension()))),
        mu((recombinationSize > 0) ? recombinationSize : lambda / 2),
        sigma0(initialStepSize) {
      }

      float_t CMAES::optimize(base::DataVector& xOpt) {
        printer.printStatusBegin("Optimizing (CMA-ES)...");

        const size_t d = f.getDimension();
        const float_t dDbl = static_cast<float_t>(d);

        base::DataVector xMean(x0);
        float_t sigma = sigma0;

        base::DataVector weights(mu);
        float_t weightsSum = 0.0;
        float_t weightsSumSquared = 0.0;

        for (size_t i = 0; i < mu; i++) {
          weights[i] = std::log((static_cast<float_t>(mu) + 0.5) /
                                static_cast<float_t>(i + 1));
          weightsSum += weights[i];
          weightsSumSquared += weights[i] * weights[i];
        }

        weights.mult(1.0 / weightsSum);
        const float_t muEff = weightsSum * weightsSum / weightsSumSquared;

        const float_t cC = (4.0 + muEff / dDbl) /
                           (dDbl + 4.0  + 2.0 * muEff / dDbl);
        const float_t cS = (muEff + 2.0) / (dDbl + muEff + 5.0);
        const float_t c1 = 2.0 / ((dDbl + 1.3) * (dDbl + 1.3) + muEff);
        const float_t cMu = std::min(
                              1.0 - c1, 2.0 * (muEff - 2.0 + 1.0 / muEff) /
                              ((dDbl + 2.0) * (dDbl + 2.0) + muEff));
        const float_t damps = 1.0 + 2.0 * std::max(0.0, std::sqrt(
                                (muEff - 1.0) / (dDbl + 1.0)) - 1.0) + cS;

        base::DataVector pC(d, 0.0);
        base::DataVector pS(d, 0.0);
        base::DataMatrix B(d, d);
        base::DataVector D(d, 1.0);
        base::DataMatrix DMatrix(d, d);

        for (size_t t1 = 0; t1 < d; t1++) {
          for (size_t t2 = 0; t2 < d; t2++) {
            B.set(t1, t2, ((t1 == t2) ? 1.0 : 0.0));
          }
        }

        base::DataMatrix C(B);
        base::DataMatrix CInvSqrt(B);
        size_t eigenEval = 0;
        const float_t chiN = std::sqrt(d) *
                             (1.0 - 1.0 / (4.0 * dDbl) +
                              1.0 / (21.0 * dDbl * dDbl));
        size_t k = 0;
        size_t countEval = 0;
        std::vector<base::DataVector> X(lambda, base::DataVector(d));
        base::DataVector fX(lambda);
        std::vector<size_t> fXIndex(lambda);
        std::vector<bool> inDomain(lambda);
        base::DataVector xMeanOld(d);

        while (countEval < N) {
          std::fill(inDomain.begin(), inDomain.end(), true);

          for (size_t t = 0; t < d; t++) {
            float_t BD_t = 0.0;

            for (size_t t2 = 0; t2 < d; t2++) {
              BD_t += B.get(t, t2) * D[t2];
            }

            for (size_t i = 0; i < lambda; i++) {
              const float_t rn = randomNumberGenerator.getGaussianRN();
              X[i][t] = xMean[t] + sigma * BD_t * rn;

              if ((X[i][t] < 0.0) || (X[i][t] > 1.0)) {
                inDomain[i] = false;
              }
            }
          }

          for (size_t i = 0; i < lambda; i++) {
            fX[i] = (inDomain[i] ? f.eval(X[i]) : INFINITY);
            fXIndex[i] = i;
          }

          countEval += lambda;
          std::sort(fXIndex.begin(), fXIndex.end(),
          [&](size_t a, size_t b) {
            return (fX.get(a) < fX.get(b));
          });
          xMeanOld = xMean;

          for (size_t t = 0; t < d; t++) {
            xMean[t] = 0.0;

            for (size_t i = 0; i < mu; i++) {
              xMean[t] += X[fXIndex[i]][t] * weights[i];
            }
          }

          for (size_t t = 0; t < d; t++) {
            float_t CInvSqrtXDiff_t = 0.0;

            for (size_t t2 = 0; t2 < d; t2++) {
              CInvSqrtXDiff_t +=
                CInvSqrt.get(t, t2) * (xMean[t2] - xMeanOld[t2]);
            }

            pS[t] = (1.0 - cS) * pS[t] + std::sqrt(cS * (2.0 - cS) * muEff) *
                    CInvSqrtXDiff_t / sigma;
          }

          const float_t pSNorm = pS.l2Norm();
          pC.mult(1.0 - cC);
          const bool hSig =
            (pSNorm / (chiN * std::sqrt(
                         1.0 - std::pow(1.0 - cS,
                                        2.0 * static_cast<float_t>(countEval) /
                                        static_cast<float_t>(lambda)))) <
             1.4 + 2.0 / (dDbl + 1.0));

          if (hSig) {
            for (size_t t = 0; t < d; t++) {
              pC[t] += std::sqrt(cC * (2.0 - cC) * muEff) *
                       (xMean[t] - xMeanOld[t]) / sigma;
            }
          }

          for (size_t t = 0; t < d; t++) {
            for (size_t t2 = 0; t2 < d; t2++) {
              const float_t CEntryOld = C.get(t, t2);
              float_t CEntryNew = (1.0 - c1 - cMu) * CEntryOld +
                                  c1 * pC[t] * pC[t2];

              if (!hSig) {
                CEntryNew += c1 * cC * (2.0 - cC) * CEntryOld;
              }

              for (size_t i = 0; i < mu; i++) {
                CEntryNew += cMu * weights[i] *
                             (X[fXIndex[i]][t] - xMeanOld[t]) *
                             (X[fXIndex[i]][t2] - xMeanOld[t2]) /
                             (sigma * sigma);
              }

              C.set(t, t2, CEntryNew);
            }
          }

          sigma *= std::exp((cS / damps) * (pSNorm / chiN - 1.0));

          if (countEval - eigenEval >
              static_cast<float_t>(lambda) / (10.0 * dDbl * (c1 + cMu))) {
            eigenEval = countEval;
            DMatrix = C;
            math::schurDecomposition(DMatrix, B);

            for (size_t t = 0; t < d; t++) {
              D[t] = std::sqrt(DMatrix.get(t, t));
            }

            for (size_t t = 0; t < d; t++) {
              for (size_t t2 = 0; t2 < d; t2++) {
                CInvSqrt.set(t, t2, B.get(t, t2) * B.get(t2, t) / D[t2]);
              }
            }
          }

          // status printing
          if (k % 10 == 0) {
            printer.printStatusUpdate(std::to_string(k) + " steps, f(x) = " +
                                      std::to_string(fX[fXIndex[0]]));
          }

          k++;

          // TODO
          /*if (D.max() > 1e7 * D.min()) {
            break;
          }*/
        }

        xOpt.resize(d);
        xOpt = X[fXIndex[0]];
        const float_t fOpt = fX[fXIndex[0]];

        printer.printStatusUpdate(std::to_string(k) + " steps, f(x) = " +
                                  std::to_string(fOpt));
        printer.printStatusEnd();

        return fOpt;
      }

      size_t CMAES::getPopulationSize() const {
        return lambda;
      }

      void CMAES::setPopulationSize(size_t populationSize) {
        lambda = populationSize;
      }

      size_t CMAES::getRecombinationSize() const {
        return mu;
      }

      void CMAES::setRecombinationSize(size_t recombinationSize) {
        mu = recombinationSize;
      }

      float_t CMAES::getInitialStepSize() const {
        return sigma0;
      }

      void CMAES::setInitialStepSize(float_t initialStepSize) {
        sigma0 = initialStepSize;
      }

    }
  }
}
