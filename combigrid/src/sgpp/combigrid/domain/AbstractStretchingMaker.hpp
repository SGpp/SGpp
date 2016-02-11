/* ****************************************************************************
 * Copyright (C) 2011 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Janos Benk (benk@in.tum.de)
#ifndef ABSTRACTSTRETCHINGMAKER_HPP_
#define ABSTRACTSTRETCHINGMAKER_HPP_

#include <vector>
#include "sgpp/combigrid/utils/combigrid_utils.hpp"

namespace combigrid {

/** Different 1D stretching's enum type */
enum Stretching { EQUIDISTANT, CHEBYSHEV, LEGENDRE, BASU, UNKNOWN, ATAN, TAN };

/**
 * enum specifying the type of coordinate transformation that has to be done to
 * map the actual coordinates to the [-1;1] interval
 *
 */
enum TRANSFORMATION_TYPE {
  FINITE,
  SEMI_INFINITE_NINF,
  SEMI_INFINITE_PINF,
  INFINITE
};

/** class to create stretching in 1D*/
class AbstractStretchingMaker {
 public:
  /**
   * @param level - integer specifying the current grid level . the
   * corresponding nr of points is 2^level + 1
   * @param min - the left boundary of the interval
   * @param max - the right boundary of the interval
   * @param stretching - the output vector of pre-computed grid points...
   * @param jacobian - the evaluated jacobian at all points of the stretching ,
   * taking into consideration
   * size of the interval and underlying tranformations.
   */

  virtual void get1DStretching(int level, double min, double max,
                               std::vector<double>* stretching,
                               std::vector<double>* jacobian) const = 0;

  virtual ~AbstractStretchingMaker() {}

  virtual Stretching getStretchingType() const = 0;

  double transforminterval(double min, double max, double point,
                           TRANSFORMATION_TYPE tp) const {
    double result = point;

    switch (tp) {
      case SEMI_INFINITE_NINF: {
        result = max - (1.0 - point) / (1.0 + point);
        break;
      }

      case SEMI_INFINITE_PINF: {
        result = min + (1.0 + point) / (1.0 - point);
        break;
      }

      case INFINITE: {
        result = point / (1.0 - point * point);
        break;
      }

      case FINITE: {
        result = 0.5 * (min + max) + 0.5 * (max - min) * point;
        break;
      }

      default:
        result = 0.5 * (min + max) + 0.5 * (max - min) * point;
        break;
    }

    return result;
  }

  double transformationJacobian(double min, double max, double pt,
                                TRANSFORMATION_TYPE type) const {
    double result = 1.0;

    // the if's are in order to avoid singularities
    switch (type) {
      case SEMI_INFINITE_NINF: {
        if (pt != -1.0) result = 2.0 / ((1.0 + pt) * (1 + pt));

        break;
      }

      case SEMI_INFINITE_PINF: {
        if (pt != 1.0) {
          result = 2 / ((1.0 - pt) * (1 - pt));
        }

        break;
      }

      case INFINITE: {
        if (pt != -1.0 && pt != 1.0) {
          result = (1 + pt * pt) / ((1.0 - pt * pt) * (1.0 - pt * pt));
        }

        break;
      }

      case FINITE: {
        result = (max - min) / 2.0;
        break;
      }

      default: {
        result = (max - min) / 2.0;
        break;
      }
    }

    return result;
  }
};
}  // namespace combigrid

#endif /* ABSTRACTSTRETCHINGMAKER_HPP_ */
