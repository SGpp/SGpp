// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef MODIFIED_POLY_BASE_HPP
#define MODIFIED_POLY_BASE_HPP

#include <cmath>
#include <vector>

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/basis/Basis.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Modified polynomial base functions.
     * Special polynomial functions to cover values unequal 0 at the border. Implemented as seen in AWR 2 paper
     * by Prof. Bungartz (http://www5.in.tum.de/wiki/index.php/Algorithmen_des_Wissenschaftlichen_Rechnens_II_-_Winter_08)
     */
    template<class LT, class IT>
    class PolyModifiedBasis: public Basis<LT, IT> {
      protected:
        /// Pointer to polynoms
        double* polynoms;
        /// polynom's max. degree
        size_t degree;

      public:

        /**
         * Constructor
         *
         * @param degree the polynom's max. degree
         */
        PolyModifiedBasis(size_t degree) : polynoms(NULL), degree(degree) {
          //if(degree < 0)
          //{
          //  throw factory_exception("PolyBasis: degree < 0");
          //}

          int polycount = (1 << (degree + 1)) - 1;
          std::vector<double> x;

          // degree + 1 for the polynom, +1 for the integral value, +1 for the base-point
          polynoms = new double[(degree + 1 + 2) * polycount];
          initPolynoms(x, 1, 1);
        }

        /**
         * Destructor
         */
        ~PolyModifiedBasis() {
          if (polynoms) {
            delete [] polynoms;
          }
        }

        /**
         * Evaluate a basis function.
         * Has a dependence on the absolute position of grid point and support.
         */
        double eval(LT level, IT index, double p) {
          size_t deg = degree + 1 < level ? degree + 1 : level;

          size_t idMask = (1 << deg) - 1;
          size_t id = (((index & idMask) >> 1) | (1 << (deg - 1))) - 1;

          // scale p to a value in [-1.0,1.0]
          double val = (1 << level) * p - index;
          return evalPolynom(id, deg, val);
        }

        double evalHierToTop(LT level, IT index, DataVector& koeffs, double pos) {
          double result = 0.0;

          for (; level >= 1; level--) {
            result += koeffs[level] * eval(level, index, pos);
            index = ((index - 1) / 2);
            index = (index % 2 == 0) ? (index + 1) : index;
          }

          return result;
        }

      private:
        /**
         * Evaluate a basis function.
         * Has a dependence on the absolute position of grid point and support.
         */
        double evalPolynom(size_t id, size_t deg, double val) {
          double* x_store = this->polynoms + (degree + 1 + 2) * id;
          double* y_store = x_store + 2;

          double y_val = y_store[deg - 1];
          // @todo (blank) scale val back into the right range
          double x_val = x_store[0] + val * pow(2.0, -(1.0) * (static_cast<double>(deg)));

          //Horner
          for (size_t i = deg - 2; i > 0; i--) {
            y_val = y_val * x_val + y_store[i];
          }

          return y_val * x_val + y_store[0];
        }

        /**
         * recursively creates polynomial values
         */
        void initPolynoms(std::vector<double>& x, LT level, IT index) {
          // Add new point
          x.push_back(index * pow(2.0, -(1.0) * (static_cast<double>(level))));

          std::vector<double> y;
          std::vector<double> intpoly;

          for (size_t i = 0; i < level - 1; i++) {
            y.push_back(0.0);
            intpoly.push_back(0.0);
          }

          y.push_back(1.0);
          intpoly.push_back(0.0);

          // Every poly has a unique id similiar to sparse grid level/index pairs
          size_t id = ((index >> 1) | (1 << (level - 1))) - 1;

          int n = level;
          std::vector<std::vector<double> > lagpoly;

          /**
           * Fill lagpoly with multiplied lagrange polynomials
           * Add lagrange polynomials together into intpoly
           */
          for (int i = 0; i < n; i++) {
            lagpoly.push_back(std::vector<double>());
            lagpoly[i].push_back(1.0);
            double fac = y[i];

            int j = 0;

            for (int k = 0; k < n; k++) {
              if (k == i) {
                continue;
              }

              lagpoly[i].push_back(lagpoly[i][j]);

              for (int jj = j; jj > 0; jj--) {
                lagpoly[i][jj] = lagpoly[i][jj - 1] - lagpoly[i][jj] * x[k];
              }

              lagpoly[i][0] *= -x[k];
              j += 1;
              fac /= (x[i] - x[k]);
            }

            for (int l = 0; l < n; l++) {
              lagpoly[i][l] *= fac;
              intpoly[l] += lagpoly[i][l];
            }
          }

          //determine position in storage. (degree + 1) polynomial factors and 2 values for integral and x-value
          double* x_store = this->polynoms + (degree + 3) * id;
          double* y_store = x_store + 2;

          // Copy values into storage
          for (int i = 0; i < n; i++) {
            y_store[i] = intpoly[i];
          }

          x_store[0] = x.back();


          if ((level) < degree + 1) {
            initPolynoms(x, level + 1, index * 2 - 1);
            initPolynoms(x, level + 1, index * 2 + 1);
          }

          x.pop_back();

        }
    };

    // default type-def (unsigned int for level and index)
    typedef PolyModifiedBasis<unsigned int, unsigned int> SPolyModifiedBase;
  }
}

#endif /* MODIFIED_POLY_BASE_HPP */