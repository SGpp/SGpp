// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/CMAES.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>

#include <algorithm>
#include <cmath>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      CMAES::CMAES(ObjectiveFunction& f, size_t maxFcnEvalCount) :
        Optimizer(f, maxFcnEvalCount) {
      }

      float_t CMAES::optimize(base::DataVector& xOpt) {
        printer.printStatusBegin("Optimizing (CMA-ES)...");

        const size_t d = f.getDimension();
        const float_t dDbl = static_cast<float_t>(d);

        base::DataVector xMean(d, 0.5); // TODO: initialize
        float_t sigma = 0.3; // TODO: initialize

        const size_t lambda = 4 + static_cast<size_t>(3.0 * std::log(d)); // TODO: initialize
        const float_t muDbl = static_cast<float_t>(lambda) / 2.0;
        const size_t mu = static_cast<size_t>(muDbl);
        base::DataVector weights(mu);
        float_t weightsSum = 0.0;
        float_t weightsSumSquared = 0.0;

        for (size_t i = 0; i < mu; i++) {
          weights[i] = std::log((muDbl + 0.5) / static_cast<float_t>(i + 1));
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
            schurDecomposition(DMatrix, B);

            for (size_t t = 0; t < d; t++) {
              D[t] = std::sqrt(DMatrix.get(t, t));
            }

            for (size_t t = 0; t < d; t++) {
              for (size_t t2 = 0; t2 < d; t2++) {
                CInvSqrt.set(t, t2, B.get(t, t2) * B.get(t2, t) / D[t2]);
              }
            }
          }

          k++;
          std::cout << "k = " << k << "\n";
          std::cout << "C = " << C.toString() << "\n";
          std::cout << "B = " << B.toString() << "\n";
          std::cout << "D = " << D.toString() << "\n";

          // TODO
          if (D.max() > 1e7 * D.min()) {
            break;
          }
        }

        xOpt.resize(d);
        xOpt = X[fXIndex[0]];
        const float_t fOpt = fX[fXIndex[0]];

        return fOpt;
      }

      void CMAES::schurDecomposition(base::DataMatrix& S,
                                     base::DataMatrix& V) {
        const size_t n = S.getNrows();
        base::DataMatrix SMinusLambdaI(n, n);
        base::DataMatrix QTilde(n, n);

        hessenbergForm(S, V);

        base::DataMatrix SNew(S);
        base::DataMatrix VNew(V);

        for (size_t k = 0; k < n - 1; k++) {
          const size_t m = n - k;
          SMinusLambdaI.resize(m, m);

          for (size_t it = 0; it < 100; it++) {
            const float_t dPlus = S.get(m - 1, m - 1) + S.get(m - 2, m - 2);
            const float_t dMinus = S.get(m - 1, m - 1) - S.get(m - 2, m - 2);
            const float_t sigma = (dMinus >= 0.0) ? 1.0 : -1.0;
            const float_t lambda = dPlus / 2.0 + sigma / 2.0 * std::sqrt(
                                     dMinus * dMinus +
                                     4.0 * S.get(m - 1, m - 2) *
                                     S.get(m - 2, m - 1));

            for (size_t i = 0; i < m; i++) {
              for (size_t j = 0; j < m; j++) {
                SMinusLambdaI.set(i, j, S.get(i, j) -
                                  ((i == j) ? lambda : 0.0));
              }
            }

            QTilde.resize(m, m);
            QRdecomposition(SMinusLambdaI, QTilde);

            for (size_t i = m; i < n; i++) {
              for (size_t j = 0; j < m; j++) {
                float_t entry = 0.0;

                for (size_t l = 0; l < m; l++) {
                  entry += S.get(i, l) * QTilde.get(l, j);
                }

                SNew.set(i, j, entry);
              }
            }

            for (size_t i = 0; i < m; i++) {
              for (size_t j = m; j < n; j++) {
                float_t entry = 0.0;

                for (size_t l = 0; l < m; l++) {
                  entry += QTilde.get(l, i) * S.get(l, j);
                }

                SNew.set(i, j, entry);
              }
            }

            for (size_t i = 0; i < m; i++) {
              for (size_t j = 0; j < m; j++) {
                float_t entry = 0.0;

                for (size_t p = 0; p < m; p++) {
                  for (size_t q = 0; q < m; q++) {
                    entry += QTilde.get(p, i) * S.get(p, q) * QTilde.get(q, j);
                  }
                }

                SNew.set(i, j, entry);
              }
            }

            S = SNew;

            for (size_t i = 0; i < n; i++) {
              for (size_t j = 0; j < m; j++) {
                float_t entry = 0.0;

                for (size_t l = 0; l < m; l++) {
                  entry += V.get(i, l) * QTilde.get(l, j);
                }

                VNew.set(i, j, entry);
              }
            }

            V = VNew;

            if (std::abs(S.get(m - 1, m - 2)) < 1e-8) {
              break;
            }
          }

          for (size_t i = m; i < n; i++) {
            S.set(i, m - 1, 0.0);
            SNew.set(i, m - 1, 0.0);
          }
        }

        for (size_t i = 1; i < n; i++) {
          S.set(i, 0, 0.0);
        }
      }

      void CMAES::hessenbergForm(base::DataMatrix& A, base::DataMatrix& V) {
        const size_t n = A.getNrows();
        base::DataMatrix Q(n, n);
        base::DataMatrix ANew(A);

        for (size_t i = 0; i < n; i++) {
          for (size_t j = 0; j < n; j++) {
            V.set(i, j, (i == j) ? 1.0 : 0.0);
          }
        }

        base::DataMatrix VNew(V);

        for (size_t k = 0; k < n - 2; k++) {
          const size_t m = n - k - 1;
          Q.resize(m, m);
          householderTransformation(A, k + 1, k, Q);

          for (size_t i = 0; i < k + 1; i++) {
            for (size_t j = k + 1; j < n; j++) {
              float_t entry = 0.0;

              for (size_t l = 0; l < m; l++) {
                entry += A.get(i, l + k + 1) * Q.get(l, j - k - 1);
              }

              ANew.set(i, j, entry);
            }
          }

          for (size_t i = k + 1; i < n; i++) {
            for (size_t j = i - 1; j < k + 1; j++) {
              float_t entry = 0.0;

              for (size_t l = 0; l < m; l++) {
                entry += Q.get(i - k - 1, l) * A.get(l + k + 1, j);
              }

              ANew.set(i, j, entry);
            }
          }

          for (size_t i = k + 2; i < n; i++) {
            ANew.set(i, k, 0.0);
          }

          for (size_t i = k + 1; i < n; i++) {
            for (size_t j = k + 1; j < n; j++) {
              float_t entry = 0.0;

              for (size_t p = 0; p < m; p++) {
                for (size_t q = 0; q < m; q++) {
                  entry += Q.get(i - k - 1, p) *
                           A.get(p + k + 1, q + k + 1) *
                           Q.get(q, j - k - 1);
                }
              }

              ANew.set(i, j, entry);
            }
          }

          A = ANew;

          for (size_t i = 0; i < n; i++) {
            for (size_t j = k + 1; j < n; j++) {
              float_t entry = 0.0;

              for (size_t l = 0; l < m; l++) {
                entry += V.get(i, l + k + 1) * Q.get(l, j - k - 1);
              }

              VNew.set(i, j, entry);
            }
          }

          V = VNew;
        }
      }

      void CMAES::QRdecomposition(base::DataMatrix& A, base::DataMatrix& Q) {
        const size_t n = A.getNrows();
        base::DataMatrix QTilde(n, n);
        base::DataVector d(n);
        base::DataMatrix ANew(A);

        for (size_t i = 0; i < n; i++) {
          for (size_t j = 0; j < n; j++) {
            Q.set(i, j, (i == j) ? 1.0 : 0.0);
          }
        }

        base::DataMatrix QNew(Q);

        for (size_t k = 0; k < n - 1; k++) {
          const size_t m = n - k;
          QTilde.resize(m, m);
          householderTransformation(A, k, k, QTilde);

          for (size_t i = k; i < n; i++) {
            for (size_t j = k; j < n; j++) {
              float_t entry = 0.0;

              for (size_t l = 0; l < m; l++) {
                entry += QTilde.get(i - k, l) * A.get(l + k, j);
              }

              ANew.set(i, j, entry);
            }
          }

          for (size_t i = k + 1; i < n; i++) {
            ANew.set(i, k, 0.0);
          }

          for (size_t i = k; i < n; i++) {
            for (size_t j = k; j < n; j++) {
              float_t entry = 0.0;

              for (size_t p = 0; p < m; p++) {
                for (size_t q = 0; q < m; q++) {
                  entry += QTilde.get(i - k, p) *
                           A.get(p + k, q + k) *
                           QTilde.get(q, j - k);
                }
              }

              ANew.set(i, j, entry);
            }
          }

          A = ANew;

          for (size_t i = 0; i < n; i++) {
            for (size_t j = k; j < n; j++) {
              float_t entry = 0.0;

              for (size_t l = 0; l < m; l++) {
                entry += Q.get(i, l + k) * QTilde.get(j - k, l);
              }

              QNew.set(i, j, entry);
            }
          }

          Q = QNew;
        }
      }

      void CMAES::householderTransformation(const base::DataMatrix& A,
                                            size_t i, size_t j,
                                            base::DataMatrix& Q) {
        const size_t n = A.getNrows();
        const size_t m = n - i;
        base::DataVector d(m);
        const float_t c1 = A.get(i, j);
        float_t cNorm = c1 * c1;

        for (size_t p = 1; p < m; p++) {
          d[p] = A.get(p + i, j);
          cNorm += d[p] * d[p];
        }

        cNorm = std::sqrt(cNorm);
        const float_t sigma = (c1 >= 0.0) ? 1.0 : -1.0;
        d[0] = c1 + sigma * cNorm;
        const float_t r = std::abs(d[0]) * cNorm;

        for (size_t p = 0; p < m; p++) {
          for (size_t q = 0; q < m; q++) {
            Q.set(p, q, ((p == q) ? 1.0 : 0.0) - d[p] * d[q] / r);
          }
        }
      }

      /*size_t DifferentialEvolution::getPopulationSize() const {
        return populationSize;
      }

      void DifferentialEvolution::setPopulationSize(size_t populationSize) {
        this->populationSize = populationSize;
      }*/

    }
  }
}
