// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/NelderMead.hpp>
#include <sgpp/optimization/tools/Permuter.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <algorithm>
#include <iostream>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      const float_t NelderMead::DEFAULT_ALPHA = 1.0;
      const float_t NelderMead::DEFAULT_BETA = 2.0;
      const float_t NelderMead::DEFAULT_GAMMA = 0.5;
      const float_t NelderMead::DEFAULT_DELTA = 0.5;
      const float_t NelderMead::STARTING_SIMPLEX_EDGE_LENGTH = 0.4;

      NelderMead::NelderMead(ObjectiveFunction& f,
                             size_t maxFcnEvalCount, float_t alpha,
                             float_t beta, float_t gamma, float_t delta) :
        Optimizer(f, maxFcnEvalCount),
        alpha(alpha),
        beta(beta),
        gamma(gamma),
        delta(delta) {
      }

      float_t NelderMead::optimize(std::vector<float_t>& xOpt) {
        printer.printStatusBegin("Optimizing (Nelder-Mead)...");

        const size_t d = f.getDimension();
        std::vector<std::vector<float_t>> points(d + 1, x0);
        std::vector<std::vector<float_t>> pointsNew(d + 1, x0);
        std::vector<float_t> fPoints(d + 1, 0.0);
        std::vector<float_t> fPointsNew(d + 1, 0.0);

        // construct starting simplex
        for (size_t t = 0; t < d; t++) {
          points[t + 1][t] = std::min(points[t + 1][t] +
                                      STARTING_SIMPLEX_EDGE_LENGTH,
                                      float_t(1.0));
          fPoints[t + 1] = f.eval(points[t + 1]);
        }

        fPoints[0] = f.eval(points[0]);

        std::vector<size_t> index(d + 1, 0);
        std::vector<float_t> pointO(d, 0.0);
        std::vector<float_t> pointR(d, 0.0);
        std::vector<float_t> pointE(d, 0.0);
        std::vector<float_t> pointIC(d, 0.0);
        std::vector<float_t> pointOC(d, 0.0);
        size_t k = 0;
        size_t numberOfFcnEvals = d + 1;

        while (true) {
          // sort points by function value
          for (size_t i = 0; i < d + 1; i++) {
            index[i] = i;
          }

          {
            Permuter<float_t> permuter(fPoints);
            std::sort(index.begin(), index.end(), permuter);
          }

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

            pointO[t] /= static_cast<float_t>(d);
            pointR[t] = pointO[t] + alpha * (pointO[t] - points[d][t]);

            if ((pointR[t] < 0.0) || (pointR[t] > 1.0)) {
              inDomain = false;
            }
          }

          float_t fPointR = (inDomain ? f.eval(pointR) : INFINITY);
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

            float_t f_point_e = (inDomain ? f.eval(pointE) : INFINITY);
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

            float_t fPointOC = (in_domain ? f.eval(pointOC) : INFINITY);
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

            float_t fPointIC = (in_domain ? f.eval(pointIC) : INFINITY);
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
                points[i][t] = points[0][t] +
                               delta * (points[i][t] - points[0][t]);

                if ((points[i][t] < 0.0) || (points[i][t] > 1.0)) {
                  in_domain = false;
                }
              }

              fPoints[i] = (in_domain ? f.eval(points[i]) : INFINITY);
            }

            numberOfFcnEvals += d;
          }

          // status printing
          if (k % 10 == 0) {
            printer.printStatusUpdate(std::to_string(k) + " steps, f(x) = " +
                                      std::to_string(fPoints[0]));
          }

          if (numberOfFcnEvals + (d + 2) > N) {
            break;
          }

          k++;
        }

        xOpt = points[0];

        printer.printStatusUpdate(std::to_string(k) + " steps, f(x) = " +
                                  std::to_string(fPoints[0]));
        printer.printStatusEnd();

        return fPoints[0];
      }

      float_t NelderMead::getAlpha() const {
        return alpha;
      }

      void NelderMead::setAlpha(float_t alpha) {
        this->alpha = alpha;
      }

      float_t NelderMead::getBeta() const {
        return beta;
      }

      void NelderMead::setBeta(float_t beta) {
        this->beta = beta;
      }

      float_t NelderMead::getGamma() const {
        return gamma;
      }

      void NelderMead::setGamma(float_t gamma) {
        this->gamma = gamma;
      }

      float_t NelderMead::getDelta() const {
        return delta;
      }

      void NelderMead::setDelta(float_t delta) {
        this->delta = delta;
      }

    }
  }
}
