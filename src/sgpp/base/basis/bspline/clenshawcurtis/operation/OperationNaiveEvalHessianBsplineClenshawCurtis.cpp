/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/basis/bspline/clenshawcurtis/operation/OperationNaiveEvalHessianBsplineClenshawCurtis.hpp"

namespace sg {
  namespace base {

    double OperationNaiveEvalHessianBsplineClenshawCurtis::evalHessian(
      DataVector& alpha, const std::vector<double>& point,
      DataVector& gradient, DataMatrix& hessian) {
      const size_t n = storage->size();
      const size_t d = storage->dim();
      double result = 0.0;

      gradient.resize(storage->dim());
      gradient.setAll(0.0);

      hessian = DataMatrix(d, d);
      hessian.setAll(0.0);

      DataVector cur_gradient(d);
      DataMatrix cur_hessian(d, d);

      for (size_t i = 0; i < n; i++) {
        const GridIndex* gp = storage->get(i);
        double cur_val = 1.0;
        cur_gradient.setAll(alpha[i]);
        cur_hessian.setAll(alpha[i]);

        for (size_t t = 0; t < d; t++) {
          double val1d = base.eval(gp->getLevel(t), gp->getIndex(t), point[t]);
          double dx1d = base.evalDx(gp->getLevel(t), gp->getIndex(t), point[t]);
          double dxdx1d = base.evalDxDx(gp->getLevel(t), gp->getIndex(t), point[t]);

          cur_val *= val1d;

          for (size_t t2 = 0; t2 < d; t2++) {
            if (t2 == t) {
              cur_gradient[t2] *= dx1d;

              for (size_t t3 = 0; t3 < d; t3++) {
                if (t3 == t) {
                  cur_hessian[t2*d+t3] *= dxdx1d;
                } else {
                  cur_hessian[t2*d+t3] *= dx1d;
                }
              }
            } else {
              cur_gradient[t2] *= val1d;

              for (size_t t3 = 0; t3 < d; t3++) {
                if (t3 == t) {
                  cur_hessian[t2*d+t3] *= dx1d;
                } else {
                  cur_hessian[t2*d+t3] *= val1d;
                }
              }
            }
          }
        }

        result += alpha[i] * cur_val;
        gradient.add(cur_gradient);
        hessian.add(cur_hessian);
      }

      return result;
    }

  }
}
