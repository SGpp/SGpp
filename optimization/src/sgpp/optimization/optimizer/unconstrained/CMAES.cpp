// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <sgpp/optimization/optimizer/unconstrained/CMAES.hpp>
#include <sgpp/optimization/tools/Math.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace sgpp {
namespace optimization {
namespace optimizer {

CMAES::CMAES(const base::ScalarFunction& f, size_t maxFcnEvalCount)
    : UnconstrainedOptimizer(f, nullptr, nullptr, maxFcnEvalCount) {}

CMAES::CMAES(const CMAES& other) : UnconstrainedOptimizer(other) {}

CMAES::~CMAES() {}

void CMAES::optimize() {
  base::Printer::getInstance().printStatusBegin("Optimizing (CMA-ES)...");

  const size_t d = f->getNumberOfParameters();

  xOpt.resize(0);
  fOpt = NAN;
  xHist.resize(0, d);
  fHist.resize(0);

  const double dDbl = static_cast<double>(d);

  const double E = std::sqrt(dDbl) * (1.0 - 1.0 / (4.0 * dDbl) + 1.0 / (21.0 * dDbl * dDbl));

  const size_t lambda = 4 + static_cast<size_t>(3.0 * std::log(dDbl));
  const double muPrime = static_cast<double>(lambda) / 2.0;
  const size_t mu = static_cast<size_t>(muPrime);

  base::DataVector wPrime(mu);

  for (size_t i = 0; i < mu; i++) {
    wPrime[i] = std::log(muPrime + 0.5) - std::log(static_cast<double>(i) + 1.0);
  }

  base::DataVector w(wPrime);
  w.mult(1.0 / wPrime.sum());

  const double muEff = 1.0 / std::pow(w.l2Norm(), 2.0);
  const double cSigma = (muEff + 2.0) / (dDbl + muEff + 5.0);
  const double dSigma =
      1.0 + 2.0 * std::max(0.0, std::sqrt((muEff - 1.0) / (dDbl + 1.0)) - 1.0) + cSigma;

  const double cC = (4.0 + muEff / dDbl) / (dDbl + 4.0 + 2.0 * muEff / dDbl);
  const double c1 = 2.0 / (std::pow(static_cast<double>(dDbl + 1.3), 2.0) + muEff);
  const double alphaMu = 2.0;
  const double cMu = std::min(1.0 - c1, alphaMu * (muEff - 2.0 + 1.0 / muEff) /
                                             (std::pow(dDbl + 2.0, 2.0) + alphaMu * muEff / 2.0));

  base::DataVector pSigma(d, 0.0);
  base::DataVector pC(d, 0.0);
  base::DataMatrix C(d, d, 0.0), CNew(d, d);
  base::DataMatrix B(d, d), D(d, d), CInvSqrt(d, d);
  base::DataVector DDiag(d);

  for (size_t t = 0; t < d; t++) {
    C(t, t) = 1.0;
  }

  base::DataVector m(x0);
  double sigma = 0.3;

  base::DataMatrix X(d, lambda), Y(d, lambda);
  base::DataVector x(d), y(d), tmp(d);
  base::DataVector fX(lambda);
  std::vector<size_t> fXOrder(lambda);

  base::DataVector yW(d);

  size_t k = 0;
  size_t numberOfFcnEvals = 0;

  while (numberOfFcnEvals < N) {
    D = C;
    math::schurDecomposition(D, B);
    D.sqrt();

    for (size_t t = 0; t < d; t++) {
      DDiag[t] = D(t, t);
    }

    for (size_t t1 = 0; t1 < d; t1++) {
      for (size_t t2 = 0; t2 < d; t2++) {
        CInvSqrt(t1, t2) = 0.0;

        for (size_t t3 = 0; t3 < d; t3++) {
          CInvSqrt(t1, t2) += B(t1, t3) * B(t2, t3) / DDiag[t3];
        }
      }
    }

    for (size_t j = 0; j < lambda; j++) {
      for (size_t t = 0; t < d; t++) {
        tmp[t] = DDiag[t] * base::RandomNumberGenerator::getInstance().getGaussianRN();
      }

      B.mult(tmp, y);
      Y.setColumn(j, y);

      x = y;
      x.mult(sigma);
      x.add(m);
      X.setColumn(j, x);

      {
        bool inDomain = true;

        for (size_t t = 0; t < d; t++) {
          if ((x[t] < 0.0) || (x[t] > 1.0)) {
            inDomain = false;
            break;
          }
        }

        fX[j] = (inDomain ? f->eval(x) : INFINITY);
        fXOrder[j] = j;
      }
    }

    numberOfFcnEvals += lambda;

    std::sort(fXOrder.begin(), fXOrder.end(),
              [&fX](size_t a, size_t b) { return (fX[a] < fX[b]); });

    for (size_t t = 0; t < d; t++) {
      yW[t] = 0.0;

      for (size_t i = 0; i < mu; i++) {
        yW[t] += w[i] * Y(t, fXOrder[i]);
      }

      m[t] += sigma * yW[t];
    }

    CInvSqrt.mult(yW, tmp);
    tmp.mult(std::sqrt(cSigma * (2.0 - cSigma) * muEff));
    pSigma.mult(1.0 - cSigma);
    pSigma.add(tmp);

    sigma *= std::exp(cSigma / dSigma * (pSigma.l2Norm() / E - 1.0));

    const double hSigma =
        ((pSigma.l2Norm() /
              std::sqrt(1.0 - std::pow(1.0 - cSigma, 2.0 * (static_cast<double>(k) + 1.0))) <
          (1.4 + 2.0 / (dDbl + 1.0)) * E)
             ? 1.0
             : 0.0);
    const double delta = (1.0 - hSigma) * cC * (2.0 - cC);

    tmp = yW;
    tmp.mult(hSigma * std::sqrt(cC * (2.0 - cC) * muEff));
    pC.mult(1.0 - cC);
    pC.add(tmp);

    for (size_t t1 = 0; t1 < d; t1++) {
      for (size_t t2 = 0; t2 < d; t2++) {
        CNew(t1, t2) = (1.0 - c1 - cMu) * C(t1, t2) + c1 * (pC[t1] * pC[t2] + delta * C(t1, t2));

        for (size_t i = 0; i < mu; i++) {
          CNew(t1, t2) += cMu * w[i] * Y(t1, fXOrder[i]) * Y(t2, fXOrder[i]);
        }
      }
    }

    C = CNew;
    k++;

    base::Printer::getInstance().printStatusUpdate(
        std::to_string(k) + " steps, f(x) = " + std::to_string(fX[fXOrder[0]]));

    X.getColumn(fXOrder[0], x);
    xHist.appendRow(x);
    fHist.append(fX[fXOrder[0]]);
  }

  xOpt.resize(d);
  X.getColumn(fXOrder[0], xOpt);
  fOpt = fX[fXOrder[0]];

  base::Printer::getInstance().printStatusUpdate(
      std::to_string(k) + " steps, f(x) = " + std::to_string(fX[fXOrder[0]]));
  base::Printer::getInstance().printStatusEnd();
}

void CMAES::clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const {
  clone = std::unique_ptr<UnconstrainedOptimizer>(new CMAES(*this));
}
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp
