// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/LaplaceEnhancedDownBBLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

LaplaceEnhancedDownBBLinear::LaplaceEnhancedDownBBLinear(sgpp::base::GridStorage* storage)
    : storage(storage),
      boundingBox(storage->getBoundingBox()),
      algoDims(storage->getAlgorithmicDimensions()),
      numAlgoDims_(storage->getAlgorithmicDimensions().size()),
      ptr_source_(nullptr),
      ptr_result_(nullptr),
      cur_algo_dim_(0),
      q_(0.0),
      t_(0.0)
#if 1
#if defined(__SSE3__)
      ,
      half_in_(_mm_set1_pd(0.5)),
      twothird_(_mm_set1_pd(2.0 / 3.0))
#endif
#else
#ifdef __SSE3__
      ,
      half_in_(_mm_set1_pd(0.5)),
      twothird_(_mm_set1_pd(2.0 / 3.0))
#endif
#endif
//      ,h_table_(nullptr),grad_table_(nullptr)
{
}

LaplaceEnhancedDownBBLinear::~LaplaceEnhancedDownBBLinear() {}

void LaplaceEnhancedDownBBLinear::operator()(sgpp::base::DataMatrix& source,
                                             sgpp::base::DataMatrix& result, grid_iterator& index,
                                             size_t dim) {
  q_ = this->boundingBox->getIntervalWidth(this->algoDims[dim]);
  t_ = this->boundingBox->getIntervalOffset(this->algoDims[dim]);

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

  if (q_ != 1.0 || t_ != 0.0) {
    size_t i = 0;

    for (i = 0; i < this->numAlgoDims_ - 1; i += 2) {
      double fl = 0.0;
      double fr = 0.0;

      if (dim == i) {
        recBB_GL(fl, fr, i, index);
      } else if (dim == i + 1) {
        recBB_LG(fl, fr, i, index);
      } else {
        recBB_LL(fl, fr, fl, fr, i, index);
      }
    }

    for (; i < this->numAlgoDims_; i++) {
      double fl = 0.0;
      double fr = 0.0;

      if (dim == i) {
        recBB_grad(i, index);
      } else {
        recBB(fl, fr, i, index);
      }
    }
  } else {
    size_t i = 0;

    for (i = 0; i < this->numAlgoDims_ - 1; i += 2) {
      double fl = 0.0;
      double fr = 0.0;
#if 1
#if defined(__SSE3__)
      __m128d fl_xmm = _mm_set1_pd(0.0);
      __m128d fr_xmm = _mm_set1_pd(0.0);
#endif
#else
#ifdef __SSE3__
      __m128d fl_xmm = _mm_set1_pd(0.0);
      __m128d fr_xmm = _mm_set1_pd(0.0);
#endif
#endif

      if (dim == i) {
        rec_GL(fl, fr, i, index);
      } else if (dim == i + 1) {
        rec_LG(fl, fr, i, index);
      } else {
#if 1
#if defined(__SSE3__)
        rec_LL(fl_xmm, fr_xmm, i, index);
#else
        rec_LL(fl, fr, fl, fr, i, index);
#endif
#else
#ifdef __SSE3__
        rec_LL(fl_xmm, fr_xmm, i, index);
#else
        rec_LL(fl, fr, fl, fr, i, index);
#endif
#endif
      }
    }

    for (; i < this->numAlgoDims_; i++) {
      double fl = 0.0;
      double fr = 0.0;

      if (dim == i) {
        rec_grad(i, index);
      } else {
        rec(fl, fr, i, index);
      }
    }
  }

  //    delete[] h_table_;
  //    delete[] grad_table_;
}

void LaplaceEnhancedDownBBLinear::rec(double fl, double fr, size_t dim, grid_iterator& index) {
  size_t seq = index.seq();
  sgpp::base::level_t l;
  sgpp::base::index_t i;
  index.get(cur_algo_dim_, l, i);

  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];

  // mesh-width
  double h = 1.0 / static_cast<double>(1 << l);

  // L2 scalar product
  double tmp_m = ((fl + fr) * 0.5);
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      (((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value)));
  double fm = tmp_m + alpha_value;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec(fl, fm, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec(fm, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::rec_grad(size_t dim, grid_iterator& index) {
  size_t seq = index.seq();
  sgpp::base::level_t l;
  sgpp::base::index_t i;
  index.get(cur_algo_dim_, l, i);

  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];

  // Gradient just in selected dimension
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      (static_cast<double>(1 << (l + 1)) * alpha_value);

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec_grad(dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec_grad(dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

#if 1
#if defined(__SSE3__)
void LaplaceEnhancedDownBBLinear::rec_LL(__m128d fl, __m128d fr, size_t dim, grid_iterator& index)
#else
void LaplaceEnhancedDownBBLinear::rec_LL(double fl, double fr, double fl2, double fr2,
                                         size_t dim, grid_iterator& index)
#endif
#else
#ifdef __SSE3__
void LaplaceEnhancedDownBBLinear::rec_LL(__m128d fl, __m128d fr, size_t dim, grid_iterator& index)
#else
void LaplaceEnhancedDownBBLinear::rec_LL(double fl, double fr, double fl2, double fr2,
                                         size_t dim, grid_iterator& index)
#endif
#endif
{
  size_t seq = index.seq();
  sgpp::base::level_t l;
  sgpp::base::index_t i;
  index.get(cur_algo_dim_, l, i);

  // mesh-width
  double h = 1.0 / static_cast<double>(1 << l);

// L2 scalar product
#if 1
#if defined(__SSE3__)
  // with intrinsics
  __m128d h_in = _mm_loaddup_pd(&h);
  __m128d fl_in = fl;
  __m128d fr_in = fr;
  __m128d alpha = _mm_loadu_pd(&ptr_source_[(seq * this->numAlgoDims_) + dim]);
  __m128d tmp = _mm_mul_pd(_mm_add_pd(fl_in, fr_in), half_in_);
  __m128d res = _mm_add_pd(_mm_mul_pd(h_in, tmp), _mm_mul_pd(alpha, _mm_mul_pd(h_in, twothird_)));
  __m128d new_fm = _mm_add_pd(alpha, tmp);
  _mm_storeu_pd(&ptr_result_[(seq * this->numAlgoDims_) + dim], res);
#else
  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  double alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];
  double tmp_m = ((fl + fr) / 2.0);
  double tmp_m2 = ((fl2 + fr2) / 2.0);
  double res = ((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value));
  double res2 = ((h * tmp_m2) + (((2.0 / 3.0) * h) * alpha_value2));
  ptr_result_[(seq * this->numAlgoDims_) + dim] = res;
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] = res2;
  double fm = tmp_m + alpha_value;
  double fm2 = tmp_m2 + alpha_value2;
#endif
#else
#ifdef __SSE3__
  // with intrinsics
  __m128d h_in = _mm_loaddup_pd(&h);
  __m128d fl_in = fl;
  __m128d fr_in = fr;
  __m128d alpha = _mm_loadu_pd(&ptr_source_[(seq * this->numAlgoDims_) + dim]);
  __m128d tmp = _mm_mul_pd(_mm_add_pd(fl_in, fr_in), half_in_);
  __m128d res = _mm_add_pd(_mm_mul_pd(h_in, tmp), _mm_mul_pd(alpha, _mm_mul_pd(h_in, twothird_)));
  __m128d new_fm = _mm_add_pd(alpha, tmp);
  _mm_storeu_pd(&ptr_result_[(seq * this->numAlgoDims_) + dim], res);
#else
  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  double alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];
  double tmp_m = ((fl + fr) / 2.0);
  double tmp_m2 = ((fl2 + fr2) / 2.0);
  double res = ((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value));
  double res2 = ((h * tmp_m2) + (((2.0 / 3.0) * h) * alpha_value2));
  ptr_result_[(seq * this->numAlgoDims_) + dim] = res;
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] = res2;
  double fm = tmp_m + alpha_value;
  double fm2 = tmp_m2 + alpha_value2;
#endif
#endif

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
#if 1
#if defined(__SSE3__)
      rec_LL(fl, new_fm, dim, index);
#else
      rec_LL(fl, fm, fl2, fm2, dim, index);
#endif
#else
#ifdef __SSE3__
      rec_LL(fl, new_fm, dim, index);
#else
      rec_LL(fl, fm, fl2, fm2, dim, index);
#endif
#endif
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
#if 1
#if defined(__SSE3__)
      rec_LL(new_fm, fr, dim, index);
#else
      rec_LL(fm, fr, fm2, fr2, dim, index);
#endif
#else
#ifdef __SSE3__
      rec_LL(new_fm, fr, dim, index);
#else
      rec_LL(fm, fr, fm2, fr2, dim, index);
#endif
#endif
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::rec_LG(double fl, double fr, size_t dim, grid_iterator& index) {
  size_t seq = index.seq();
  sgpp::base::level_t l;
  sgpp::base::index_t i;
  index.get(cur_algo_dim_, l, i);

  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  double alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];

  // mesh-width
  double h = 1.0 / static_cast<double>(1 << l);

  // L2 scalar product
  double tmp_m = ((fl + fr) * 0.5);
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      (((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value)));
  // Gradient in second dimension
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] =
      (static_cast<double>(1 << (l + 1)) * alpha_value2);

  double fm = tmp_m + alpha_value;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec_LG(fl, fm, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec_LG(fm, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::rec_GL(double fl, double fr, size_t dim, grid_iterator& index) {
  size_t seq = index.seq();
  sgpp::base::level_t l;
  sgpp::base::index_t i;
  index.get(cur_algo_dim_, l, i);

  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  double alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];

  // mesh-width
  double h = 1.0 / static_cast<double>(1 << l);

  double tmp_m = ((fl + fr) * 0.5);

  // Gradient in second dimension
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      ((static_cast<double>(1 << (l + 1))) * alpha_value);
  // L2 scalar product
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] =
      (((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value2)));

  double fm = tmp_m + alpha_value2;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec_GL(fl, fm, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec_GL(fm, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::recBB(double fl, double fr, size_t dim, grid_iterator& index) {
  size_t seq = index.seq();
  sgpp::base::level_t l;
  sgpp::base::index_t i;
  index.get(cur_algo_dim_, l, i);

  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];

  // mesh-width
  double h = 1.0 / static_cast<double>(1 << l);

  // L2 scalar product
  double tmp_m = ((fl + fr) * 0.5);
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      (((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value)) * q_);
  double fm = tmp_m + alpha_value;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB(fl, fm, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB(fm, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::recBB_LL(double fl, double fr, double fl2, double fr2,
                                           size_t dim, grid_iterator& index) {
  size_t seq = index.seq();
  sgpp::base::level_t l;
  sgpp::base::index_t i;
  index.get(cur_algo_dim_, l, i);

  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  double alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];

  // mesh-width
  double h = 1.0 / static_cast<double>(1 << l);

  // L2 scalar product
  double tmp_m = ((fl + fr) * 0.5);
  double tmp_m2 = ((fl2 + fr2) * 0.5);
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      (((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value)) * q_);
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] =
      (((h * tmp_m2) + (((2.0 / 3.0) * h) * alpha_value2)) * q_);
  double fm = tmp_m + alpha_value;
  double fm2 = tmp_m2 + alpha_value2;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB_LL(fl, fm, fl2, fm2, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB_LL(fm, fr, fm2, fr2, dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::recBB_LG(double fl, double fr, size_t dim,
                                           grid_iterator& index) {
  size_t seq = index.seq();
  sgpp::base::level_t l;
  sgpp::base::index_t i;
  index.get(cur_algo_dim_, l, i);

  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  double alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];

  // mesh-width
  double h = 1.0 / static_cast<double>(1 << l);

  // L2 scalar product
  double tmp_m = ((fl + fr) * 0.5);
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      (((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value)) * q_);
  // Gradient in second dimension
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] =
      ((static_cast<double>(1 << (l + 1)) / q_) * alpha_value2);

  double fm = tmp_m + alpha_value;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB_LG(fl, fm, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB_LG(fm, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::recBB_GL(double fl, double fr, size_t dim,
                                           grid_iterator& index) {
  size_t seq = index.seq();
  sgpp::base::level_t l;
  sgpp::base::index_t i;
  index.get(cur_algo_dim_, l, i);

  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  double alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];

  // mesh-width
  double h = 1.0 / static_cast<double>(1 << l);

  double tmp_m = ((fl + fr) * 0.5);
  // Gradient in second dimension
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      ((static_cast<double>(1 << (l + 1)) / q_) * alpha_value);
  // L2 scalar product
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] =
      (((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value2)) * q_);
  double fm = tmp_m + alpha_value2;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB_GL(fl, fm, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB_GL(fm, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

void LaplaceEnhancedDownBBLinear::recBB_grad(size_t dim, grid_iterator& index) {
  size_t seq = index.seq();
  sgpp::base::level_t l;
  sgpp::base::index_t i;
  index.get(cur_algo_dim_, l, i);

  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];

  // Gradient just in selected dimension
  ptr_result_[(seq * this->numAlgoDims_) + dim] =
      ((static_cast<double>(1 << (l + 1)) / q_) * alpha_value);

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB_grad(dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB_grad(dim, index);
    }

    index.up(cur_algo_dim_);
  }
}

}  // namespace pde
}  // namespace sgpp
