/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_SLE_SYSTEM_SYSTEM_HPP
#define SGPP_OPT_SLE_SYSTEM_SYSTEM_HPP

#include <vector>
#include <cstddef>

namespace sg {
  namespace opt {
    namespace sle {
      namespace system {

        /**
         * Abstract class representing a system of linear equations.
         * All row and column indices are zero based.
         */
        class System {
          public:
            /**
             * Constructor.
             */
            System() {
            }

            /**
             * Virtual destructor.
             */
            virtual ~System() {
            }

            /**
             * Pure virtual method for checking if a matrix entry vanishes or not.
             *
             * @param i     row index
             * @param j     column index
             * @return      whether the (i,j)-th entry of the matrix is non-zero
             */
            virtual bool isMatrixEntryNonZero(size_t i, size_t j) = 0;

            /**
             * Pure virtual method for retrieving a matrix entry.
             *
             * @param i     row index
             * @param j     column index
             * @return      (i,j)-th entry of the matrix
             */
            virtual double getMatrixEntry(size_t i, size_t j) = 0;

            /**
             * Multiply the matrix with a vector.
             * Standard implementation with \f$\mathcal{O}(n^2)\f$ scalar multiplications.
             *
             * @param       x   vector to be multiplied
             * @param[out]  y   \f$y = Ax\f$
             */
            virtual void matrixVectorMultiplication(const std::vector<double>& x,
                                                    std::vector<double>& y) {
              const size_t n = getDimension();
              y = std::vector<double>(n, 0.0);

              for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                  y[i] += getMatrixEntry(i, j) * x[j];
                }
              }
            }

            /**
             * Count all non-zero entries.
             * Standard implementation with \f$\mathcal{O}(n^2)\f$ checks.
             */
            virtual size_t countNNZ() {
              const size_t n = getDimension();
              size_t nnz = 0;

              for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                  if (isMatrixEntryNonZero(i, j)) {
                    nnz++;
                  }
                }
              }

              return nnz;
            }

            /**
             * Pure virtual method returning the dimension (number of rows/columns) of the system.
             *
             * @return  system dimension
             */
            virtual size_t getDimension() const = 0;

            /**
             * @return whether this system derives from Cloneable or not (standard: false)
             */
            virtual bool isCloneable() const {
              return false;
            }
        };

      }
    }
  }
}

#endif
