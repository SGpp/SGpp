/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/grid/type/LinearGrid.hpp"
#include "base/grid/generation/GridGenerator.hpp"

#include "parallel/pde/basis/linear/noboundary/operation/OperationLaplaceVectorizedLinear.hpp"

#include <cmath>

namespace sg {
  namespace parallel {

    OperationLaplaceVectorizedLinear::OperationLaplaceVectorizedLinear(sg::base::GridStorage* storage) : storage(storage) {
      this->level_ = new sg::base::DataMatrix(storage->size(), storage->dim());
      this->level_int_ = new sg::base::DataMatrix(storage->size(), storage->dim());
      this->index_ = new sg::base::DataMatrix(storage->size(), storage->dim());
      lcl_q = new double[this->storage->dim()];
      lcl_q_inv = new double[this->storage->dim()];

      storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
      storage->getLevelForIntegral(*(this->level_int_));
      this->lambda_ = new sg::base::DataVector(storage->dim());
      this->lambda_->setAll(1.0);
    }

    OperationLaplaceVectorizedLinear::OperationLaplaceVectorizedLinear(sg::base::GridStorage* storage, sg::base::DataVector& lambda) : storage(storage) {
      this->level_ = new sg::base::DataMatrix(storage->size(), storage->dim());
      this->level_int_ = new sg::base::DataMatrix(storage->size(), storage->dim());
      this->index_ = new sg::base::DataMatrix(storage->size(), storage->dim());
      lcl_q = new double[this->storage->dim()];
      lcl_q_inv = new double[this->storage->dim()];

      storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
      storage->getLevelForIntegral(*(this->level_int_));
      this->lambda_ = new sg::base::DataVector(lambda);
    }

    OperationLaplaceVectorizedLinear::~OperationLaplaceVectorizedLinear() {
      delete this->level_;
      delete this->level_int_;
      delete this->index_;
      delete[] lcl_q;
      delete[] lcl_q_inv;

      if (this->lambda_ != NULL)
        delete this->lambda_;
    }

    double OperationLaplaceVectorizedLinear::gradient(size_t i, size_t j, size_t dim) {
      double grad;

      double i_level_grad = level_->get(i, dim);
      double i_index_grad = index_->get(i, dim);
      double j_level_grad = level_->get(j, dim);
      double j_index_grad = index_->get(j, dim);

      //only affects the diagonal of the stiffness matrix
      bool doGrad = ((i_level_grad == j_level_grad) && (i_index_grad == j_index_grad));
      grad = i_level_grad * 2.0 * static_cast<double>(static_cast<int>(doGrad));

      return (grad * lcl_q_inv[dim]);
    }

    double OperationLaplaceVectorizedLinear::l2dot(size_t i, size_t j, size_t dim) {
      double lid = level_->get(i, dim);
      double ljd = level_->get(j, dim);
      double iid = index_->get(i, dim);
      double ijd = index_->get(j, dim);
      double in_lid = level_int_->get(i, dim);
      double in_ljd = level_int_->get(j, dim);

      //////////
      /// First case: lid == ljd
      /// we mask this in the end in the last line
      //////////

      // use textbook formular if both operands are identical
      // ansatz function on the same level but with different indecies
      // don't overlap!
      double res_one = (2.0 / 3.0) * in_lid * (iid == ijd);

      //////////
      /// Second case: lid != ljd
      /// we mask this in the end in the last line
      //////////

      // now we select the 1st as the "narrow" basisfunction (it has a higher level)
      // --> we know can regard 2nd function as linear function and therefore
      // apply the wellknown formular: 1/2 * (f_l + f_r) * 2^{-l}
      bool selector = (lid > ljd);
      double i1d = iid * (selector) + ijd * (!selector);
      //double l1d = lid*(selector) + ljd*(!selector);
      double in_l1d = in_lid * (selector) + in_ljd * (!selector);
      double i2d = ijd * (selector) + iid * (!selector);
      double l2d = ljd * (selector) + lid * (!selector);
      double in_l2d = in_ljd * (selector) + in_lid * (!selector);

      // check if Ansatz functions on different
      // levels do not overlap and neg.
      // overlap is 1 if the functions overlap
      // if they don't overlap result is zero.
      // we mask the l2 scalar product in the end
      double q = (i1d - 1) * in_l1d;
      double p = (i1d + 1) * in_l1d;
      bool overlap = (std::max(q, (i2d - 1) * in_l2d) < std::min(p, (i2d + 1) * in_l2d));

#define ALEX
#ifdef BENNI
      // We determine the distance and mirrow the
      // it to the left, the we apply the gradient of the
      // linear overlapped part. Finally we transform
      // back the our actual basis function by
      // multiplying with 2^{-l}
      double diff = (i1d * in_l1d) - (i2d * in_l2d);
      double temp_res = fabs(diff - in_l1d) + fabs(diff + in_l1d) - fabs(diff);
      temp_res *= l2d;
      temp_res = (1 - temp_res) * in_l1d;
#endif
#ifdef ALEX
      // we determine fl and fr by plugging them
      // into the sparse grids basis functions given by l2d and i2d.
      // Then we use the formular from above: 1/2 * (f_l + f_r) * 2^{-l}
      double temp_res = 2.0 - fabs(l2d * q - i2d) - fabs(l2d * p - i2d);
      temp_res *= (0.5 * in_l1d);
#endif
      double res_two = temp_res * overlap; // Now mask result

      return (res_one * (lid == ljd) + res_two * (lid != ljd)) * lcl_q[dim];
    }

    void OperationLaplaceVectorizedLinear::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      result.setAll(0.0);

      // fill q array
      for (size_t d = 0; d < this->storage->dim(); d++) {
        sg::base::BoundingBox* boundingBox = this->storage->getBoundingBox();
        lcl_q[d] = boundingBox->getIntervalWidth(d);
        lcl_q_inv[d] = 1.0 / boundingBox->getIntervalWidth(d);
      }

#if 0
      #pragma omp parallel
      {

        double* gradient_temp = new double[this->storage->dim()];
        double* dot_temp = new double[this->storage->dim()];

        #pragma omp for

        for (size_t i = 0; i < this->storage->size(); i++) {
          for (size_t j = 0; j < this->storage->size(); j++) {
            for (size_t d = 0; d < this->storage->dim(); d++) {
              gradient_temp[d] = gradient(i, j, d);
            }

            for (size_t d = 0; d < this->storage->dim(); d++) {
              dot_temp[d] = l2dot(i, j, d);
            }

            for (size_t d_outer = 0; d_outer < this->storage->dim(); d_outer++) {
              double element = alpha[j];

              for (size_t d_inner = 0; d_inner < this->storage->dim(); d_inner++) {
                //element *= ((l2dot(i, j, d_inner)*(d_outer != d_inner)) + (gradient(i, j, d_inner)*(d_outer == d_inner)));
                element *= ((dot_temp[d_inner] * (d_outer != d_inner)) + (gradient_temp[d_inner] * (d_outer == d_inner)));
              }

              result[i] += element;
            }
          }
        }

        delete [] gradient_temp;
        delete [] dot_temp;

      }
#elif 1

      //        size_t gbl_count=0;

      #pragma omp parallel
      {
        double* gradient_temp = new double[this->storage->dim()];
        double* dot_temp = new double[this->storage->dim()];
        double all_gradient_zero = 0;
        //  size_t lcl_count = 0;

        #pragma omp for

        for (size_t ii = 0; ii < this->storage->size(); ii++) {
          for (size_t jj = 0; jj < this->storage->size(); jj++) {
            all_gradient_zero = 0;

            for (size_t d = 0; d < this->storage->dim(); d++) {
              gradient_temp[d] = gradient(ii, jj, d);
              all_gradient_zero += gradient_temp[d];
            }

            if (all_gradient_zero > 0) {
              for (size_t d = 0; d < this->storage->dim(); d++) {
                dot_temp[d] = l2dot(ii, jj, d);
              }

              for (size_t d_outer = 0; d_outer < this->storage->dim(); d_outer++) {
                double element = alpha[jj];

                for (size_t d_inner = 0; d_inner < this->storage->dim(); d_inner++) {
                  //element *= ((l2dot(i, j, d_inner)*(d_outer != d_inner)) + (gradient(i, j, d_inner)*(d_outer == d_inner)));
                  element *= ((dot_temp[d_inner] * (d_outer != d_inner)) + (gradient_temp[d_inner] * (d_outer == d_inner)));
                }

                result[ii] += (this->lambda_->get(d_outer) * element);
              }

              //        lcl_count++;
            }
          }
        }

        delete [] gradient_temp;
        delete [] dot_temp;

        //  #pragma omp critical
        //  {
        //    gbl_count += lcl_count;
        //  }
        //
      }
      //  size_t no_basisfunctions = (this->storage->size()*this->storage->size());
      //  std::cout << "skipped evaluations: " <<  no_basisfunctions - gbl_count << " savings: " << ((double)(no_basisfunctions - gbl_count))/((double)no_basisfunctions)*100 << "%" << std::endl;
#else
      //Original Version
      #pragma omp parallel for

      for (size_t i = 0; i < this->storage->size(); i++) {
        for (size_t j = 0; j < this->storage->size(); j++) {
          for (size_t d_outer = 0; d_outer < this->storage->dim(); d_outer++) {
            double element = alpha[j];

            for (size_t d_inner = 0; d_inner < this->storage->dim(); d_inner++) {
              //element *= ((l2dot(i, j, d_inner)*(d_outer != d_inner)) + (gradient(i, j, d_inner)*(d_outer == d_inner)));
              element *= ((l2dot(i, j, d_inner) * (d_outer != d_inner)) + (gradient(i, j, d_inner) * (d_outer == d_inner)));
            }

            result[i] += (this->lambda_->get(d_outer) * element);
          }
        }
      }

#endif
    }

  }

}
