// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TOOLS_MATH_HPP
#define SGPP_OPTIMIZATION_TOOLS_MATH_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Namespace with linear algebra functions.
     */
    namespace math {

      /**
       * Schur decomposition of given matrix \f$A\f$ with transformation
       * matrix \f$V\f$
       * (similiarity transformation such that \f$S = V^{-1} AV\f$
       * upper triangular, \f$V\f$ orthogonal).
       * The Schur decomposition exists iff the characteristic polynomial
       * of \f$A\f$ factorizes in real linear factors.
       *
       * @param[in,out] square matrix, afterwards \f$S\f$
       * @param[out]    transformation matrix
       *                (must have correct size \f$n \times n\f$)
       */
      void schurDecomposition(base::DataMatrix& A, base::DataMatrix& V);

      /**
       * QR decomposition of given matrix \f$A\f$ with transformation
       * matrix \f$Q\f$ with \f$A = QR\f$
       * (\f$Q\f$ orthogonal and \f$R\f$ upper triangular).
       *
       * @param[in,out] square matrix, afterwards \f$R\f$
       * @param[out]    transformation matrix
       *                (must have correct size \f$n \times n\f$)
       */
      void QRDecomposition(base::DataMatrix& A, base::DataMatrix& Q);

      /**
       * Hessenberg form of given matrix \f$A\f$ with transformation
       * matrix \f$V\f$ with \f$H = V^{-1} AV\f$ in Hessenberg form
       * (similiarity transformation such that entries \f$(i,j)\f$ vanish
       * for \f$i > j + 1\f$, \f$V\f$ orthogonal).
       *
       * @param[in,out] square matrix, afterwards \f$H\f$
       * @param[out]    transformation matrix
       *                (must have correct size \f$n \times n\f$)
       */
      void hessenbergForm(base::DataMatrix& A, base::DataMatrix& V);

      /**
       * Calculate transformation matrix \f$Q\f$ of a Householder
       * transformation.
       * The normal vector (defining the reflection hyperplane) used is
       * \f$d := (c_1 + \sigma \lVert c \rVert_2, c_2, \dotsc, c_n)\f$
       * with \f$c := A(i:end,j)\f$ and \f$\sigma := +1\f$ for
       * \f$\sigma \ge 0\f$ and \f$\sigma := -1\f$ otherwise.
       * \f$Q\f$ is symmetric and orthogonal
       * (i.e. \f$Q = Q^{\mathrm{t}} = Q^{-1}\f$).
       * After applying \f$QA\f$, the entries \f$(i+1:end,j)\f$
       * should vanish.
       *
       * @param A       square matrix containing the normal vector
       * @param i       row index of starting row of normal vector
       * @param j       column index of starting column of normal vector
       * @param[out] Q  transformation matrix
       *                (must have correct size \f$(n-i) \times (n-i)\f$)
       */
      void householderTransformation(const base::DataMatrix& A,
                                     size_t i, size_t j,
                                     base::DataMatrix& Q);

    }
  }
}

#endif /* SGPP_OPTIMIZATION_TOOLS_MATH_HPP */
