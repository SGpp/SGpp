/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_SLE_SYSTEM_FULL_HPP
#define SGPP_OPT_SLE_SYSTEM_FULL_HPP

#include "opt/sle/system/Cloneable.hpp"
#include "base/datatypes/DataMatrix.hpp"

#include <vector>
#include <cstddef>

namespace sg {
  namespace opt {
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
            inline double getMatrixEntry(size_t i, size_t j) {
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
             * @return smart pointer to cloned object
             */
            virtual tools::SmartPointer<Cloneable> clone() {
              return tools::SmartPointer<Cloneable>(new Full(A));
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
