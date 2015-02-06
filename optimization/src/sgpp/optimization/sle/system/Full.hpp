// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_SLE_SYSTEM_FULL_HPP
#define SGPP_OPTIMIZATION_SLE_SYSTEM_FULL_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/system/Cloneable.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <vector>
#include <cstddef>

namespace SGPP {
  namespace optimization {
    namespace sle {
      namespace system {

        /**
         * Full linear system, essentially a wrapper around base::DataMatrix.
         */
        class Full : public Cloneable {
          public:
            /**
             * Constructor.
             * Do not destruct the matrix A before this object!
             *
             * @param A     coefficient matrix
             */
            Full(base::DataMatrix& A) : Cloneable(), A(A) {
            }

            /**
             * Virtual destructor.
             */
            virtual ~Full() {
            }

            /**
             * @param i     row index
             * @param j     column index
             * @return      whether the (i,j)-th entry of the matrix is non-zero
             */
            inline bool isMatrixEntryNonZero(size_t i, size_t j) {
              return (A.get(i, j) != 0.0);
            }

            /**
             * @param i     row index
             * @param j     column index
             * @return      (i,j)-th entry of the matrix
             */
            inline float_t getMatrixEntry(size_t i, size_t j) {
              return A.get(i, j);
            }

            /**
             * @return  coefficient matrix
             */
            base::DataMatrix& getA() {
              return A;
            }

            /**
             * @param A coefficient matrix (do not destruct before this object!)
             */
            void setA(base::DataMatrix& A) {
              this->A = A;
            }

            size_t getDimension() const {
              return A.getNrows();
            }

            /**
             * Clones the linear system.
             * Because A is stored as a reference, A is not copied (only b).
             *
             * @param[out] clone pointer to cloned object
             */
            virtual void clone(Cloneable*& clone) {
              clone = new Full(A);
            }

          protected:
            /// coefficient matrix
            base::DataMatrix& A;
        };

      }
    }
  }
}

#endif
