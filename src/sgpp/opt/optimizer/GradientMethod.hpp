/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_OPTIMIZER_GRADIENTMETHOD_HPP
#define SGPP_OPT_OPTIMIZER_GRADIENTMETHOD_HPP

#include "opt/optimizer/Optimizer.hpp"
#include "opt/function/ObjectiveGradient.hpp"

namespace sg {
  namespace opt {
    namespace optimizer {

      /**
       * Gradient-based method of steepest descent.
       */
      class GradientMethod : public Optimizer {
        public:
          /// default maximal number of iterations
          static const size_t DEFAULT_MAX_IT_COUNT = 2000;
          /// default beta (parameter for Armijo's rule)
          static const double DEFAULT_BETA;
          /// default gamma (parameter for Armijo's rule)
          static const double DEFAULT_GAMMA;
          /// default tolerance (parameter for Armijo's rule)
          static const double DEFAULT_TOLERANCE;
          /// default epsilon (parameter for Armijo's rule)
          static const double DEFAULT_EPSILON;

          /**
           * Constructor.
           *
           * @param f             objective function
           * @param f_gradient    objective function gradient
           * @param max_it_count  maximal number of iterations
           * @param beta          beta (parameter for Armijo's rule)
           * @param gamma         gamma (parameter for Armijo's rule)
           * @param tolerance     tolerance (parameter for Armijo's rule)
           * @param epsilon       epsilon (parameter for Armijo's rule)
           */
          GradientMethod(function::Objective& f,
                         function::ObjectiveGradient& f_gradient,
                         size_t max_it_count = DEFAULT_MAX_IT_COUNT,
                         double beta = DEFAULT_BETA,
                         double gamma = DEFAULT_GAMMA,
                         double tolerance = DEFAULT_TOLERANCE,
                         double epsilon = DEFAULT_EPSILON);

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
           * @return smart pointer to objective function gradient
           */
          const tools::SmartPointer<function::ObjectiveGradient>& getObjectiveGradient() const;

          /**
           * @return              beta (parameter for Armijo's rule)
           */
          double getBeta() const;

          /**
           * @param beta          beta (parameter for Armijo's rule)
           */
          void setBeta(double beta);

          /**
           * @return              gamma (parameter for Armijo's rule)
           */
          double getGamma() const;

          /**
           * @param gamma         gamma (parameter for Armijo's rule)
           */
          void setGamma(double gamma);

          /**
           * @return              tolerance (parameter for Armijo's rule)
           */
          double getTolerance() const;

          /**
           * @param tolerance     tolerance (parameter for Armijo's rule)
           */
          void setTolerance(double tolerance);

          /**
           * @return              epsilon (parameter for Armijo's rule)
           */
          double getEpsilon() const;

          /**
           * @param epsilon       epsilon (parameter for Armijo's rule)
           */
          void setEpsilon(double epsilon);

        protected:
          /// objective function gradient
          tools::SmartPointer<function::ObjectiveGradient> f_gradient;
          /// beta (parameter for Armijo's rule)
          double beta;
          /// gamma (parameter for Armijo's rule)
          double gamma;
          /// tolerance (parameter for Armijo's rule)
          double tol;
          /// epsilon (parameter for Armijo's rule)
          double eps;
      };

    }
  }
}

#endif
