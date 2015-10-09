// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Math.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
  namespace optimization {
    namespace math {

      void schurDecomposition(base::DataMatrix& A, base::DataMatrix& V) {
        const size_t n = A.getNrows();
        base::DataMatrix AMinusLambdaI(n, n);
        base::DataMatrix QTilde(n, n);

        // calculate Hessenberg form of A
        hessenbergForm(A, V);

        base::DataMatrix ANew(A);
        base::DataMatrix VNew(V);

        // for each eigenvalue
        for (size_t k = 0; k < n - 1; k++) {
          // find eigenvalue of A(1:m,1:m)
          const size_t m = n - k;
          AMinusLambdaI.resize(m, m);

          // do QR iteration
          for (size_t it = 0; it < 100; it++) {
            // break loop if lowest subdiagonal entry is small enough
            // ==> convergence
            if (std::abs(A(m - 1, m - 2)) < 1e-8) {
              break;
            }

            // calculate shift: lambda is the eigenvalue of A(m-1:m,m-1:m)
            // which is closest to A(m,m)
            const float_t dPlus = A(m - 1, m - 1) + A(m - 2, m - 2);
            const float_t dMinus = A(m - 1, m - 1) - A(m - 2, m - 2);
            const float_t sigma = (dMinus >= 0.0) ? 1.0 : -1.0;
            const float_t lambda = dPlus / 2.0 + sigma / 2.0 * std::sqrt(
                                     dMinus * dMinus +
                                     4.0 * A(m - 1, m - 2) * A(m - 2, m - 1));

            // calculate A(1:m,1:m) - lambda*eye(m)
            for (size_t i = 0; i < m; i++) {
              for (size_t j = 0; j < m; j++) {
                AMinusLambdaI(i, j) = A(i, j) - ((i == j) ? lambda : 0.0);
              }
            }

            // calculate QR decomposition of AMinusLambdaI
            QTilde.resize(m, m);
            QRDecomposition(AMinusLambdaI, QTilde);

            // we want now to update
            // A = (A_ul A_ur) --> (Qt 0)^T * A * (Qt 0) = (Qt^T*A_ul*Qt  Qt^T*A_ur)
            //     (A_bl A_br)     (0  I)         (0  I)   (     A_bl*Qt       A_br)
            // and
            // V = (V_l V_r) --> V * (Qt 0) = (V_l*Qt  V_r)
            //                       (0  I)
            // with I being the identity matrix and Qt = QTilde

            // update lower left values
            for (size_t i = m; i < n; i++) {
              for (size_t j = 0; j < m; j++) {
                float_t entry = 0.0;

                for (size_t l = 0; l < m; l++) {
                  entry += A(i, l) * QTilde(l, j);
                }

                ANew(i, j) = entry;
              }
            }

            // update upper right values
            for (size_t i = 0; i < m; i++) {
              for (size_t j = m; j < n; j++) {
                float_t entry = 0.0;

                for (size_t l = 0; l < m; l++) {
                  entry += QTilde(l, i) * A(l, j);
                }

                ANew(i, j) = entry;
              }
            }

            // update upper left values
            for (size_t i = 0; i < m; i++) {
              for (size_t j = 0; j < m; j++) {
                float_t entry = 0.0;

                for (size_t p = 0; p < m; p++) {
                  for (size_t q = 0; q < m; q++) {
                    entry += QTilde(p, i) * A(p, q) * QTilde(q, j);
                  }
                }

                ANew(i, j) = entry;
              }
            }

            A = ANew;

            // update transformation matrix (only left part)
            for (size_t i = 0; i < n; i++) {
              for (size_t j = 0; j < m; j++) {
                float_t entry = 0.0;

                for (size_t l = 0; l < m; l++) {
                  entry += V(i, l) * QTilde(l, j);
                }

                VNew(i, j) = entry;
              }
            }

            V = VNew;
          }

          // set lower entries explicitly to zero
          for (size_t i = m; i < n; i++) {
            A(i, m - 1) = 0.0;
            ANew(i, m - 1) = 0.0;
          }
        }

        // set lower entries explicitly to zero
        for (size_t i = 1; i < n; i++) {
          A(i, 0) = 0.0;
        }
      }

      void QRDecomposition(base::DataMatrix& A, base::DataMatrix& Q) {
        const size_t n = A.getNrows();
        base::DataMatrix QTilde(n, n);
        base::DataVector d(n);
        base::DataMatrix ANew(A);

        // set transformation matrix to the identity matrix
        for (size_t i = 0; i < n; i++) {
          for (size_t j = 0; j < n; j++) {
            Q(i, j) = ((i == j) ? 1.0 : 0.0);
          }
        }

        base::DataMatrix QNew(Q);

        // generate zeros in the k-th column
        for (size_t k = 0; k < n - 1; k++) {
          // size of current Householder transformation
          const size_t m = n - k;
          QTilde.resize(m, m);
          // calculate Householder transformation
          householderTransformation(A, k, k, QTilde);

          // we want now to update
          // A = (A_ul A_ur) --> (I  0) * A = (A_ul     A_ur)
          //     (   0 A_br)     (0 Qt)       (   0  Qt*A_br)
          // and
          // Q = (Q_l Q_r) --> Q * (I  0)^T = (Q_l  Q_r*Qt^T)
          //                       (0 Qt)
          // with I being the identity matrix and Qt = QTilde

          // update lower right values (lower left part is zero)
          for (size_t i = k; i < n; i++) {
            for (size_t j = k; j < n; j++) {
              float_t entry = 0.0;

              for (size_t l = 0; l < m; l++) {
                entry += QTilde(i - k, l) * A(l + k, j);
              }

              ANew(i, j) = entry;
            }
          }

          // set values exactly to zero
          // (guaranteed by Householder transformation)
          for (size_t i = k + 1; i < n; i++) {
            ANew(i, k) = 0.0;
          }

          A = ANew;

          // update transformation matrix (only right part)
          for (size_t i = 0; i < n; i++) {
            for (size_t j = k; j < n; j++) {
              float_t entry = 0.0;

              for (size_t l = 0; l < m; l++) {
                entry += Q(i, l + k) * QTilde(j - k, l);
              }

              QNew(i, j) = entry;
            }
          }

          Q = QNew;
        }
      }

      void hessenbergForm(base::DataMatrix& A, base::DataMatrix& V) {
        const size_t n = A.getNrows();
        base::DataMatrix QTilde(n, n);
        base::DataMatrix ANew(A);

        // set transformation matrix to the identity matrix
        for (size_t i = 0; i < n; i++) {
          for (size_t j = 0; j < n; j++) {
            V(i, j) = ((i == j) ? 1.0 : 0.0);
          }
        }

        base::DataMatrix VNew(V);

        // generate zeros in the k-th column
        for (size_t k = 0; k < n - 2; k++) {
          // size of current Householder transformation
          const size_t m = n - k - 1;
          QTilde.resize(m, m);
          // calculate Householder transformation
          householderTransformation(A, k + 1, k, QTilde);

          // we want now to update
          // A = (A_ul A_ur) --> (I  0) * A * (I  0) = (   A_ul     A_ur*Qt)
          //     (A_bl A_br)     (0 Qt)       (0 Qt)   (Qt*A_bl  Qt*A_br*Qt)
          // and
          // V = (V_l V_r) --> V * (I  0) = (V_l  V_r*Qt)
          //                       (0 Qt)
          // with I being the identity matrix and Qt = QTilde

          // update upper right values
          for (size_t i = 0; i < k + 1; i++) {
            for (size_t j = k + 1; j < n; j++) {
              float_t entry = 0.0;

              for (size_t l = 0; l < m; l++) {
                entry += A(i, l + k + 1) * QTilde(l, j - k - 1);
              }

              ANew(i, j) = entry;
            }
          }

          // update lower left values
          for (size_t i = k + 1; i < n; i++) {
            for (size_t j = i - 1; j < k + 1; j++) {
              float_t entry = 0.0;

              for (size_t l = 0; l < m; l++) {
                entry += QTilde(i - k - 1, l) * A(l + k + 1, j);
              }

              ANew(i, j) = entry;
            }
          }

          // set values exactly to zero
          // (guaranteed by Householder transformation)
          for (size_t i = k + 2; i < n; i++) {
            ANew(i, k) = 0.0;
          }

          // update lower right values
          for (size_t i = k + 1; i < n; i++) {
            for (size_t j = k + 1; j < n; j++) {
              float_t entry = 0.0;

              for (size_t p = 0; p < m; p++) {
                for (size_t q = 0; q < m; q++) {
                  entry += QTilde(i - k - 1, p) *
                           A(p + k + 1, q + k + 1) *
                           QTilde(q, j - k - 1);
                }
              }

              ANew(i, j) = entry;
            }
          }

          A = ANew;

          // update transformation matrix (only right part)
          for (size_t i = 0; i < n; i++) {
            for (size_t j = k + 1; j < n; j++) {
              float_t entry = 0.0;

              for (size_t l = 0; l < m; l++) {
                entry += V(i, l + k + 1) * QTilde(l, j - k - 1);
              }

              VNew(i, j) = entry;
            }
          }

          V = VNew;
        }
      }

      void householderTransformation(const base::DataMatrix& A,
                                     size_t i, size_t j,
                                     base::DataMatrix& Q) {
        // dimension of matrix
        const size_t n = A.getNrows();
        // dimension of transformation
        const size_t m = n - i;
        // normal vector to reflection hyperplane
        base::DataVector d(m);

        const float_t c1 = A(i, j);
        float_t cNorm = c1 * c1;

        // read from matrix
        for (size_t p = 1; p < m; p++) {
          d[p] = A(p + i, j);
          cNorm += d[p] * d[p];
        }

        cNorm = std::sqrt(cNorm);
        const float_t sigma = (c1 >= 0.0) ? 1.0 : -1.0;
        d[0] = c1 + sigma * cNorm;
        const float_t r = std::abs(d[0]) * cNorm;

        // create transformation matrix
        for (size_t p = 0; p < m; p++) {
          for (size_t q = 0; q < m; q++) {
            Q(p, q) = ((p == q) ? 1.0 : 0.0) - d[p] * d[q] / r;
          }
        }
      }

    }
  }
}
