// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/LaplaceEnhancedDownBBLinearBoundary.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

LaplaceEnhancedDownBBLinearBoundary::LaplaceEnhancedDownBBLinearBoundary(
    sgpp::base::GridStorage* storage)
    : LaplaceEnhancedDownBBLinear(storage) {}

LaplaceEnhancedDownBBLinearBoundary::~LaplaceEnhancedDownBBLinearBoundary() {}

void LaplaceEnhancedDownBBLinearBoundary::operator()(sgpp::base::DataMatrix& source,
                                                     sgpp::base::DataMatrix& result,
                                                     grid_iterator& index, size_t dim) {
  q_ = this->boundingBox->getIntervalWidth(this->algoDims[dim]);
  t_ = this->boundingBox->getIntervalOffset(this->algoDims[dim]);
  double q_reci = 1.0 / q_;

  //    h_table_ = new double[MAX_TABLE_DEPTH+1];
  //    grad_table_ = new double[MAX_TABLE_DEPTH+1];
  //
  //    for (int i = 0; i <= MAX_TABLE_DEPTH; i++)
  //    {
  //        h_table_[i] = 1.0/static_cast<double>(1<<i);
  //        grad_table_[i] = (static_cast<double>(1<<(i+1)))/q;
  //    }

  ptr_source_ = source.getPointer();
  ptr_result_ = result.getPointer();
  cur_algo_dim_ = this->algoDims[dim];

  /*
   * Handle Level 0
   */
  // This handles the diagonal only
  //////////////////////////////////////
  // left boundary
  index.resetToLeftLevelZero(this->algoDims[dim]);
  size_t seq_left = index.seq();

  // right boundary
  index.resetToRightLevelZero(this->algoDims[dim]);
  size_t seq_right = index.seq();

  if (q_ != 1.0 || t_ != 0.0) {
    size_t i = 0;

    for (i = 0; i < this->numAlgoDims_ - 1; i += 2) {
      double fl1 = source.get(seq_left, i);
      double fl2 = source.get(seq_left, i + 1);
      double fr1 = source.get(seq_right, i);
      double fr2 = source.get(seq_right, i + 1);

      if (dim == i) {
        calcGradBoundary(fl1, fr1, seq_left, seq_right, i, cur_algo_dim_, q_reci);
        calcL2Boundary(fl2, fr2, seq_left, seq_right, i + 1, cur_algo_dim_, q_);

        if (!index.hint()) {
          index.resetToLevelOne(dim);

          if (!this->storage->isInvalidSequenceNumber(index.seq())) {
            recBB_GL(fl2, fr2, i, index);
          }

          index.resetToLeftLevelZero(dim);
        }
      } else if (dim == i + 1) {
        calcL2Boundary(fl1, fr1, seq_left, seq_right, i, cur_algo_dim_, q_);
        calcGradBoundary(fl2, fr2, seq_left, seq_right, i + 1, cur_algo_dim_, q_reci);

        if (!index.hint()) {
          index.resetToLevelOne(dim);

          if (!this->storage->isInvalidSequenceNumber(index.seq())) {
            recBB_LG(fl1, fr1, i, index);
          }

          index.resetToLeftLevelZero(dim);
        }
      } else {
        calcL2Boundary(fl1, fr1, seq_left, seq_right, i, cur_algo_dim_, q_);
        calcL2Boundary(fl2, fr2, seq_left, seq_right, i + 1, cur_algo_dim_, q_);

        if (!index.hint()) {
          index.resetToLevelOne(dim);

          if (!this->storage->isInvalidSequenceNumber(index.seq())) {
            recBB_LL(fl1, fr1, fl2, fr2, i, index);
          }

          index.resetToLeftLevelZero(dim);
        }
      }
    }

    for (; i < this->numAlgoDims_; i++) {
      double fl = source.get(seq_left, i);
      double fr = source.get(seq_right, i);

      if (dim == i) {
        calcGradBoundary(fl, fr, seq_left, seq_right, i, cur_algo_dim_, q_reci);

        if (!index.hint()) {
          index.resetToLevelOne(dim);

          if (!this->storage->isInvalidSequenceNumber(index.seq())) {
            recBB_grad(i, index);
          }

          index.resetToLeftLevelZero(dim);
        }
      } else {
        calcL2Boundary(fl, fr, seq_left, seq_right, i, cur_algo_dim_, q_);

        if (!index.hint()) {
          index.resetToLevelOne(dim);

          if (!this->storage->isInvalidSequenceNumber(index.seq())) {
            recBB(fl, fr, i, index);
          }

          index.resetToLeftLevelZero(dim);
        }
      }
    }
  } else {
    size_t i = 0;

    for (i = 0; i < this->numAlgoDims_ - 1; i += 2) {
      double fl1 = source.get(seq_left, i);
      double fl2 = source.get(seq_left, i + 1);
      double fr1 = source.get(seq_right, i);
      double fr2 = source.get(seq_right, i + 1);
#if 1
#if defined(__SSE3__)
      __m128d fl_xmm = _mm_set_pd(fl2, fl1);
      __m128d fr_xmm = _mm_set_pd(fr2, fr1);
#endif
#else
#ifdef __SSE3__
      __m128d fl_xmm = _mm_set_pd(fl2, fl1);
      __m128d fr_xmm = _mm_set_pd(fr2, fr1);
#endif
#endif

      if (dim == i) {
        calcGradBoundary(fl1, fr1, seq_left, seq_right, i, cur_algo_dim_, 1.0);
        calcL2Boundary(fl2, fr2, seq_left, seq_right, i + 1, cur_algo_dim_, 1.0);

        if (!index.hint()) {
          index.resetToLevelOne(dim);

          if (!this->storage->isInvalidSequenceNumber(index.seq())) {
            rec_GL(fl2, fr2, i, index);
          }

          index.resetToLeftLevelZero(dim);
        }
      } else if (dim == i + 1) {
        calcL2Boundary(fl1, fr1, seq_left, seq_right, i, cur_algo_dim_, 1.0);
        calcGradBoundary(fl2, fr2, seq_left, seq_right, i + 1, cur_algo_dim_, 1.0);

        if (!index.hint()) {
          index.resetToLevelOne(dim);

          if (!this->storage->isInvalidSequenceNumber(index.seq())) {
            rec_LG(fl1, fr1, i, index);
          }

          index.resetToLeftLevelZero(dim);
        }
      } else {
        calcL2Boundary(fl1, fr1, seq_left, seq_right, i, cur_algo_dim_, 1.0);
        calcL2Boundary(fl2, fr2, seq_left, seq_right, i + 1, cur_algo_dim_, 1.0);
#if 1
#if defined(__SSE3__)

        if (!index.hint()) {
          index.resetToLevelOne(dim);

          if (!this->storage->isInvalidSequenceNumber(index.seq())) {
            rec_LL(fl_xmm, fr_xmm, i, index);
          }

          index.resetToLeftLevelZero(dim);
        }

#else

        if (!index.hint()) {
          index.resetToLevelOne(dim);

          if (!this->storage->isInvalidSequenceNumber(index.seq())) {
            rec_LL(fl1, fr1, fl2, fr2, i, index);
          }

          index.resetToLeftLevelZero(dim);
        }

#endif
#else
#ifdef __SSE3__

        if (!index.hint()) {
          index.resetToLevelOne(dim);

          if (!this->storage->isInvalidSequenceNumber(index.seq())) {
            rec_LL(fl_xmm, fr_xmm, i, index);
          }

          index.resetToLeftLevelZero(dim);
        }

#else

        if (!index.hint()) {
          index.resetToLevelOne(dim);

          if (!this->storage->isInvalidSequenceNumber(index.seq())) {
            rec_LL(fl1, fr1, fl2, fr2, i, index);
          }

          index.resetToLeftLevelZero(dim);
        }

#endif
#endif
      }
    }

    for (; i < this->numAlgoDims_; i++) {
      double fl = source.get(seq_left, i);
      double fr = source.get(seq_right, i);

      if (dim == i) {
        calcGradBoundary(fl, fr, seq_left, seq_right, i, cur_algo_dim_, 1.0);

        if (!index.hint()) {
          index.resetToLevelOne(dim);

          if (!this->storage->isInvalidSequenceNumber(index.seq())) {
            rec_grad(i, index);
          }

          index.resetToLeftLevelZero(dim);
        }
      } else {
        calcL2Boundary(fl, fr, seq_left, seq_right, i, cur_algo_dim_, 1.0);

        if (!index.hint()) {
          index.resetToLevelOne(dim);

          if (!this->storage->isInvalidSequenceNumber(index.seq())) {
            rec(fl, fr, i, index);
          }

          index.resetToLeftLevelZero(dim);
        }
      }
    }
  }
}

void LaplaceEnhancedDownBBLinearBoundary::calcL2Boundary(double fl, double fr, size_t seq_left,
                                                         size_t seq_right, size_t dim,
                                                         size_t algo_dim, double q) {
  if (this->boundingBox->hasDirichletBoundaryLeft(algo_dim))
    ptr_result_[(seq_left * this->numAlgoDims_) + dim] = 0.0;
  else
    ptr_result_[(seq_left * this->numAlgoDims_) + dim] = (1.0 / 3.0) * fl * q;

  if (this->boundingBox->hasDirichletBoundaryRight(algo_dim))
    ptr_result_[(seq_right * this->numAlgoDims_) + dim] = 0.0;
  else
    ptr_result_[(seq_right * this->numAlgoDims_) + dim] =
        ((1.0 / 3.0) * fr * q) + ((1.0 / 6.0) * fl * q);
}

void LaplaceEnhancedDownBBLinearBoundary::calcGradBoundary(double fl, double fr, size_t seq_left,
                                                           size_t seq_right, size_t dim,
                                                           size_t algo_dim, double q_reci) {
  if (this->boundingBox->hasDirichletBoundaryLeft(algo_dim))
    ptr_result_[(seq_left * this->numAlgoDims_) + dim] = 0.0;
  else
    ptr_result_[(seq_left * this->numAlgoDims_) + dim] = fl;

  if (this->boundingBox->hasDirichletBoundaryRight(cur_algo_dim_))
    ptr_result_[(seq_right * this->numAlgoDims_) + dim] = 0.0;
  else
    ptr_result_[(seq_right * this->numAlgoDims_) + dim] = fr + ((-1.0 * q_reci) * fl);
}

}  // namespace pde
}  // namespace sgpp
