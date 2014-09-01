/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_OPTIMIZER_NELDERMEAD_HPP
#define SGPP_OPT_OPTIMIZER_NELDERMEAD_HPP

#include "opt/optimizer/Optimizer.hpp"

namespace sg {
  namespace opt {
    namespace optimizer {

      /**
       * Gradient-free Nelder-Mead method.
       */
      class NelderMead : public Optimizer {
        public:
          /// default reflection coefficient
          static const double DEFAULT_ALPHA;
          /// default expansion coefficient
          static const double DEFAULT_BETA;
          /// default contraction coefficient
          static const double DEFAULT_GAMMA;
          /// default shrinking coefficient
          static const double DEFAULT_DELTA;
          /// default maximal number of function evaluations
          static const size_t DEFAULT_MAX_FCN_EVAL_COUNT = 1000;
          /// edge length of starting simplex
          static const double STARTING_SIMPLEX_EDGE_LENGTH;

          /**
           * Constructor.
           * The starting point is set to \f$(0.5, \dotsc, 0.5)^{\mathrm{T}}\f$.
           *
           * @param f                     objective function
           * @param max_fcn_eval_count    maximal number of function evaluations
           * @param alpha                 reflection coefficient
           * @param beta                  expansion coefficient
           * @param gamma                 contraction coefficient
           * @param delta                 shrinking coefficient
           */
          NelderMead(function::Objective& f,
                     size_t max_fcn_eval_count = DEFAULT_MAX_FCN_EVAL_COUNT,
                     double alpha = DEFAULT_ALPHA,
                     double beta = DEFAULT_BETA,
                     double gamma = DEFAULT_GAMMA,
                     double delta = DEFAULT_DELTA);

          /**
           * @param[out] xopt optimal point
           * @return          optimal objective function value
           */
          double optimize(std::vector<double>& xopt);

          /**
           * @return smart pointer to cloned object
           */
          tools::SmartPointer<Optimizer> clone();

          /**
           * @return          reflection coefficient
           */
          double getAlpha() const;

          /**
           * @param alpha     reflection coefficient
           */
          void setAlpha(double alpha);

          /**
           * @return          expansion coefficient
           */
          double getBeta() const;

          /**
           * @param beta      expansion coefficient
           */
          void setBeta(double beta);

          /**
           * @return          contraction coefficient
           */
          double getGamma() const;

          /**
           * @param gamma     contraction coefficient
           */
          void setGamma(double gamma);

          /**
           * @return          shrinking coefficient
           */
          double getDelta() const;

          /**
           * @param delta     shrinking coefficient
           */
          void setDelta(double delta);

        protected:
          /// reflection coefficient
          double alpha;
          /// expansion coefficient
          double beta;
          /// contraction coefficient
          double gamma;
          /// shrinking coefficient
          double delta;
      };

    }
  }
}

#endif
