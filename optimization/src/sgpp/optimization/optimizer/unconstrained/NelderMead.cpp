// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/optimization/optimizer/unconstrained/NelderMead.hpp>

#include <algorithm>
#include <limits>
#include <vector>

namespace sgpp {
namespace optimization {
namespace optimizer {

NelderMead::NelderMead(const base::ScalarFunction& f, size_t maxFcnEvalCount, double alpha,
                       double beta, double gamma, double delta)
    : UnconstrainedOptimizer(f, nullptr, nullptr, maxFcnEvalCount),
      alpha(alpha),
      beta(beta),
      gamma(gamma),
      delta(delta) {}

NelderMead::NelderMead(const NelderMead& other)
    : UnconstrainedOptimizer(other),
      alpha(other.alpha),
      beta(other.beta),
      gamma(other.gamma),
      delta(other.delta) {}

NelderMead::~NelderMead() {}

void NelderMead::optimize() {
  base::Printer::getInstance().printStatusBegin("Optimizing (Nelder-Mead)...");

  const size_t d = f->getNumberOfParameters();

  xOpt.resize(0);
  fOpt = NAN;
  xHist.resize(0, d);
  fHist.resize(0);

  std::vector<base::DataVector> points(d + 1, x0);
  std::vector<base::DataVector> pointsNew(d + 1, x0);
  base::DataVector fPoints(d + 1);
  base::DataVector fPointsNew(d + 1);

  // construct starting simplex
  for (size_t t = 0; t < d; t++) {
    points[t + 1][t] = std::min(points[t + 1][t] + STARTING_SIMPLEX_EDGE_LENGTH, 1.0);
    fPoints[t + 1] = f->eval(points[t + 1]);
  }

  fPoints[0] = f->eval(points[0]);

  std::vector<size_t> index(d + 1, 0);
  base::DataVector pointO(d);
  base::DataVector pointR(d);
  base::DataVector pointE(d);
  base::DataVector pointIC(d);
  base::DataVector pointOC(d);
  size_t k = 0;
  size_t numberOfFcnEvals = d + 1;

  while (true) {
    // sort points by function value
    for (size_t i = 0; i < d + 1; i++) {
      index[i] = i;
    }

    std::sort(index.begin(), index.end(),
              [&fPoints](size_t a, size_t b) { return (fPoints[a] < fPoints[b]); });

    // that could be solved more efficiently, but it suffices for now
    for (size_t i = 0; i < d + 1; i++) {
      pointsNew[i] = points[index[i]];
      fPointsNew[i] = fPoints[index[i]];
    }

    points = pointsNew;
    fPoints = fPointsNew;

    bool inDomain = true;
    bool shrink = false;

    // calculate point_o (barycenter of all points but the last) and
    // point_r (reflected point) simultaneously
    for (size_t t = 0; t < d; t++) {
      pointO[t] = 0.0;

      for (size_t i = 0; i < d; i++) {
        pointO[t] += points[i][t];
      }

      pointO[t] /= static_cast<double>(d);
      pointR[t] = pointO[t] + alpha * (pointO[t] - points[d][t]);

      if ((pointR[t] < 0.0) || (pointR[t] > 1.0)) {
        inDomain = false;
      }
    }

    double fPointR = (inDomain ? f->eval(pointR) : std::numeric_limits<double>::infinity());
    numberOfFcnEvals++;

    if ((fPoints[0] <= fPointR) && (fPointR < fPoints[d - 1])) {
      points[d] = pointR;
      fPoints[d] = fPointR;
    } else if (fPointR < fPoints[0]) {
      bool inDomain = true;

      // calculate expanded point
      for (size_t t = 0; t < d; t++) {
        pointE[t] = pointO[t] + beta * (pointR[t] - pointO[t]);

        if ((pointE[t] < 0.0) || (pointE[t] > 1.0)) {
          inDomain = false;
        }
      }

      double f_point_e = (inDomain ? f->eval(pointE) : std::numeric_limits<double>::infinity());
      numberOfFcnEvals++;

      if (f_point_e < fPointR) {
        points[d] = pointE;
        fPoints[d] = f_point_e;
      } else {
        points[d] = pointR;
        fPoints[d] = fPointR;
      }
    } else if (fPointR < fPoints[d]) {
      bool in_domain = true;

      // calculate outer contracted point
      for (size_t t = 0; t < d; t++) {
        pointOC[t] = pointO[t] + gamma * (pointR[t] - pointO[t]);

        if ((pointOC[t] < 0.0) || (pointOC[t] > 1.0)) {
          in_domain = false;
        }
      }

      double fPointOC = (in_domain ? f->eval(pointOC) : std::numeric_limits<double>::infinity());
      numberOfFcnEvals++;

      if (fPointOC <= fPointR) {
        points[d] = pointOC;
        fPoints[d] = fPointOC;
      } else {
        shrink = true;
      }
    } else {
      bool in_domain = true;

      // calculate inner contracted point
      for (size_t t = 0; t < d; t++) {
        pointIC[t] = pointO[t] - gamma * (pointO[t] - points[d][t]);

        if ((pointIC[t] < 0.0) || (pointIC[t] > 1.0)) {
          in_domain = false;
        }
      }

      double fPointIC = (in_domain ? f->eval(pointIC) : std::numeric_limits<double>::infinity());
      numberOfFcnEvals++;

      if (fPointIC < fPoints[d]) {
        points[d] = pointIC;
        fPoints[d] = fPointIC;
      } else {
        shrink = true;
      }
    }

    if (shrink) {
      // shrink all points but the first
      for (size_t i = 1; i < d + 1; i++) {
        bool in_domain = true;

        for (size_t t = 0; t < d; t++) {
          points[i][t] = points[0][t] + delta * (points[i][t] - points[0][t]);

          if ((points[i][t] < 0.0) || (points[i][t] > 1.0)) {
            in_domain = false;
          }
        }

        fPoints[i] = (in_domain ? f->eval(points[i]) : std::numeric_limits<double>::infinity());
      }

      numberOfFcnEvals += d;
    }

    // status printing
    if (k % 10 == 0) {
      base::Printer::getInstance().printStatusUpdate(
          std::to_string(k) + " steps, f(x) = " + std::to_string(fPoints[0]));
    }

    if (numberOfFcnEvals + (d + 2) > N) {
      break;
    }

    k++;
    xHist.appendRow(points[0]);
    fHist.append(fPoints[0]);
  }

  xOpt.resize(d);
  xOpt = points[0];
  fOpt = fPoints[0];

  base::Printer::getInstance().printStatusUpdate(std::to_string(k) +
                                                 " steps, f(x) = " + std::to_string(fPoints[0]));
  base::Printer::getInstance().printStatusEnd();
}

double NelderMead::getAlpha() const { return alpha; }

void NelderMead::setAlpha(double alpha) { this->alpha = alpha; }

double NelderMead::getBeta() const { return beta; }

void NelderMead::setBeta(double beta) { this->beta = beta; }

double NelderMead::getGamma() const { return gamma; }

void NelderMead::setGamma(double gamma) { this->gamma = gamma; }

double NelderMead::getDelta() const { return delta; }

void NelderMead::setDelta(double delta) { this->delta = delta; }

void NelderMead::clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const {
  clone = std::unique_ptr<UnconstrainedOptimizer>(new NelderMead(*this));
}
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp
