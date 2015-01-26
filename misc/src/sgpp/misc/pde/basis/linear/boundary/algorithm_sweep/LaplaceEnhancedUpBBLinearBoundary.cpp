/* *****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                          *
* This file is part of the SG++ project. For conditions of distribution and    *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp           *
***************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/misc/pde/basis/linear/boundary/algorithm_sweep/LaplaceEnhancedUpBBLinearBoundary.hpp>

namespace sg {
  namespace pde {

    LaplaceEnhancedUpBBLinearBoundary::LaplaceEnhancedUpBBLinearBoundary(sg::base::GridStorage* storage) : LaplaceEnhancedUpBBLinear(storage) {
    }

    LaplaceEnhancedUpBBLinearBoundary::~LaplaceEnhancedUpBBLinearBoundary() {
    }

    void LaplaceEnhancedUpBBLinearBoundary::operator()(sg::base::DataMatrix& source, sg::base::DataMatrix& result, grid_iterator& index, size_t dim) {
      q_ = this->boundingBox->getIntervalWidth(this->algoDims[dim]);
      t_ = this->boundingBox->getIntervalOffset(this->algoDims[dim]);
      double q_reci = 1.0 / q_;

      ptr_source_ = source.getPointer();
      ptr_result_ = result.getPointer();
      cur_algo_dim_ = this->algoDims[dim];

      size_t seq_left;
      size_t seq_right;

      // left boundary
      seq_left = index.seq();
      // right boundary
      index.right_levelzero(this->algoDims[dim]);
      seq_right = index.seq();
      index.left_levelzero(this->algoDims[dim]);

      if (q_ != 1.0 || t_ != 0.0) {
        size_t i = 0;

        for (i = 0; i < this->numAlgoDims_ - 1; i += 2) {
          if (dim == i) {
            double fl2 = 0.0;
            double fr2 = 0.0;

            if (!index.hint()) {
              index.top(dim);

              if (!this->storage->end(index.seq())) {
                recBB_GL(fl2, fr2, i, index);
              }

              index.left_levelzero(dim);
            }

            calcGradBoundary(0.0, 0.0, seq_left, seq_right, i, cur_algo_dim_, q_reci);
            calcL2Boundary(fl2, fr2, seq_left, seq_right, i + 1, cur_algo_dim_, q_);
          } else if (dim == i + 1) {
            double fl1 = 0.0;
            double fr1 = 0.0;

            if (!index.hint()) {
              index.top(dim);

              if (!this->storage->end(index.seq())) {
                recBB_LG(fl1, fr1, i, index);
              }

              index.left_levelzero(dim);
            }

            calcL2Boundary(fl1, fr1, seq_left, seq_right, i, cur_algo_dim_, q_);
            calcGradBoundary(0.0, 0.0, seq_left, seq_right, i + 1, cur_algo_dim_, q_reci);
          } else {
            double fl1 = 0.0;
            double fr1 = 0.0;
            double fl2 = 0.0;
            double fr2 = 0.0;

            if (!index.hint()) {
              index.top(dim);

              if (!this->storage->end(index.seq())) {
                recBB_LL(fl1, fr1, fl2, fr2, i, index);
              }

              index.left_levelzero(dim);
            }

            calcL2Boundary(fl1, fr1, seq_left, seq_right, i, cur_algo_dim_, q_);
            calcL2Boundary(fl2, fr2, seq_left, seq_right, i + 1, cur_algo_dim_, q_);
          }
        }

        for ( ; i < this->numAlgoDims_; i++) {
          if (dim == i) {
            if (!index.hint()) {
              index.top(dim);

              if (!this->storage->end(index.seq())) {
                recBB_grad(i, index);
              }

              index.left_levelzero(dim);
            }

            calcGradBoundary(0.0, 0.0, seq_left, seq_right, i, cur_algo_dim_, q_reci);
          } else {
            double fl = 0.0;
            double fr = 0.0;

            if (!index.hint()) {
              index.top(dim);

              if (!this->storage->end(index.seq())) {
                recBB(fl, fr, i, index);
              }

              index.left_levelzero(dim);
            }

            calcL2Boundary(fl, fr, seq_left, seq_right, i, cur_algo_dim_, q_);
          }
        }
      } else {
        size_t i = 0;

        for (i = 0; i < this->numAlgoDims_ - 1; i += 2) {
          if (dim == i) {
            double fl2 = 0.0;
            double fr2 = 0.0;

            if (!index.hint()) {
              index.top(dim);

              if (!this->storage->end(index.seq())) {
                rec_GL(fl2, fr2, i, index);
              }

              index.left_levelzero(dim);
            }

            calcGradBoundary(0.0, 0.0, seq_left, seq_right, i, cur_algo_dim_, 1.0);
            calcL2Boundary(fl2, fr2, seq_left, seq_right, i + 1, cur_algo_dim_, 1.0);
          } else if (dim == i + 1) {
            double fl1 = 0.0;
            double fr1 = 0.0;

            if (!index.hint()) {
              index.top(dim);

              if (!this->storage->end(index.seq())) {
                rec_LG(fl1, fr1, i, index);
              }

              index.left_levelzero(dim);
            }

            calcL2Boundary(fl1, fr1, seq_left, seq_right, i, cur_algo_dim_, 1.0);
            calcGradBoundary(0.0, 0.0, seq_left, seq_right, i + 1, cur_algo_dim_, 1.0);
          } else {
            double fl1 = 0.0;
            double fr1 = 0.0;
            double fl2 = 0.0;
            double fr2 = 0.0;
#ifdef __SSE3__
            __m128d fl_xmm = _mm_set1_pd(0.0);
            __m128d fr_xmm = _mm_set1_pd(0.0);

            if (!index.hint()) {
              index.top(dim);

              if (!this->storage->end(index.seq())) {
                rec_LL(fl_xmm, fr_xmm, i, index);
              }

              index.left_levelzero(dim);
            }

            _mm_storel_pd(&fl1, fl_xmm);
            _mm_storeh_pd(&fl2, fl_xmm);
            _mm_storel_pd(&fr1, fr_xmm);
            _mm_storeh_pd(&fr2, fr_xmm);
#else

            if (!index.hint()) {
              index.top(dim);

              if (!this->storage->end(index.seq())) {
                rec_LL(fl1, fr1, fl2, fr2, i, index);
              }

              index.left_levelzero(dim);
            }

#endif
            calcL2Boundary(fl1, fr1, seq_left, seq_right, i, cur_algo_dim_, 1.0);
            calcL2Boundary(fl2, fr2, seq_left, seq_right, i + 1, cur_algo_dim_, 1.0);
          }
        }

        for ( ; i < this->numAlgoDims_; i++) {
          if (dim == i) {
            if (!index.hint()) {
              index.top(dim);

              if (!this->storage->end(index.seq())) {
                rec_grad(i, index);
              }

              index.left_levelzero(dim);
            }

            calcGradBoundary(0.0, 0.0, seq_left, seq_right, i, cur_algo_dim_, 1.0);
          } else {
            double fl = 0.0;
            double fr = 0.0;

            if (!index.hint()) {
              index.top(dim);

              if (!this->storage->end(index.seq())) {
                rec(fl, fr, i, index);
              }

              index.left_levelzero(dim);
            }

            calcL2Boundary(fl, fr, seq_left, seq_right, i, cur_algo_dim_, 1.0);
          }
        }
      }
    }

    void LaplaceEnhancedUpBBLinearBoundary::calcL2Boundary(double fl, double fr, size_t seq_left, size_t seq_right, size_t dim, size_t algo_dim, double q) {
      if (this->boundingBox->hasDirichletBoundaryLeft(algo_dim))
        ptr_result_[(seq_left * this->numAlgoDims_) + dim] = 0.0;
      else
        ptr_result_[(seq_left * this->numAlgoDims_) + dim] =  fl + ((1.0 / 6.0) * ptr_source_[(seq_right * this->numAlgoDims_) + dim] * q);

      if (this->boundingBox->hasDirichletBoundaryRight(algo_dim))
        ptr_result_[(seq_right * this->numAlgoDims_) + dim] = 0.0;
      else
        ptr_result_[(seq_right * this->numAlgoDims_) + dim] = fr;
    }

    void LaplaceEnhancedUpBBLinearBoundary::calcGradBoundary(double fl, double fr, size_t seq_left, size_t seq_right, size_t dim, size_t algo_dim, double q_reci) {
      if (this->boundingBox->hasDirichletBoundaryLeft(algo_dim))
        ptr_result_[(seq_left * this->numAlgoDims_) + dim] = 0.0;
      else
        ptr_result_[(seq_left * this->numAlgoDims_) + dim] = ((-1.0 * q_reci) * ptr_source_[(seq_right * this->numAlgoDims_) + dim]);

      if (this->boundingBox->hasDirichletBoundaryRight(cur_algo_dim_))
        ptr_result_[(seq_right * this->numAlgoDims_) + dim] = 0.0;
      else
        ptr_result_[(seq_right * this->numAlgoDims_) + dim] = 0.0;
    }

    // namespace pde
  }
  // namespace sg
}
