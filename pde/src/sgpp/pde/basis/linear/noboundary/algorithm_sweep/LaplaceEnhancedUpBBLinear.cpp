// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/LaplaceEnhancedUpBBLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

LaplaceEnhancedUpBBLinear::LaplaceEnhancedUpBBLinear(sgpp::base::GridStorage* storage)
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
      half_in_(_mm_set1_pd(0.5))
#endif
#else
#ifdef __SSE3__
      ,
      half_in_(_mm_set1_pd(0.5))
#endif
#endif
{
}

LaplaceEnhancedUpBBLinear::~LaplaceEnhancedUpBBLinear() {}

void LaplaceEnhancedUpBBLinear::operator()(sgpp::base::DataMatrix& source,
                                           sgpp::base::DataMatrix& result, grid_iterator& index,
                                           size_t dim) {
  q_ = this->boundingBox->getIntervalWidth(this->algoDims[dim]);
  t_ = this->boundingBox->getIntervalOffset(this->algoDims[dim]);

  ptr_source_ = source.getPointer();
  ptr_result_ = result.getPointer();
  cur_algo_dim_ = this->algoDims[dim];

  if (q_ != 1.0 || t_ != 0.0) {
    size_t i = 0;

    for (i = 0; i < this->numAlgoDims_ - 1; i += 2) {
      if (dim == i) {
        double fl = 0.0;
        double fr = 0.0;
        recBB_GL(fl, fr, i, index);
      } else if (dim == i + 1) {
        double fl = 0.0;
        double fr = 0.0;
        recBB_LG(fl, fr, i, index);
      } else {
        double fl = 0.0;
        double fr = 0.0;
        double fl2 = 0.0;
        double fr2 = 0.0;
        recBB_LL(fl, fr, fl2, fr2, i, index);
      }
    }

    for (; i < this->numAlgoDims_; i++) {
      if (dim == i) {
        recBB_grad(i, index);
      } else {
        double fl = 0.0;
        double fr = 0.0;
        recBB(fl, fr, i, index);
      }
    }
  } else {
    size_t i = 0;

    for (i = 0; i < this->numAlgoDims_ - 1; i += 2) {
      if (dim == i) {
        double fl = 0.0;
        double fr = 0.0;
        rec_GL(fl, fr, i, index);
      } else if (dim == i + 1) {
        double fl = 0.0;
        double fr = 0.0;
        rec_LG(fl, fr, i, index);
      } else {
#if 1
#if defined(__SSE3__)
        __m128d fl_xmm = _mm_set1_pd(0.0);
        __m128d fr_xmm = _mm_set1_pd(0.0);

        rec_LL(fl_xmm, fr_xmm, i, index);
#else
        double fl = 0.0;
        double fr = 0.0;
        double fl2 = 0.0;
        double fr2 = 0.0;

        rec_LL(fl, fr, fl2, fr2, i, index);
#endif
#else
#ifdef __SSE3__
        __m128d fl_xmm = _mm_set1_pd(0.0);
        __m128d fr_xmm = _mm_set1_pd(0.0);

        rec_LL(fl_xmm, fr_xmm, i, index);
#else
        double fl = 0.0;
        double fr = 0.0;
        double fl2 = 0.0;
        double fr2 = 0.0;

        rec_LL(fl, fr, fl2, fr2, i, index);
#endif
#endif
      }
    }

    for (; i < this->numAlgoDims_; i++) {
      if (dim == i) {
        rec_grad(i, index);
      } else {
        double fl = 0.0;
        double fr = 0.0;
        rec(fl, fr, i, index);
      }
    }
  }
}

void LaplaceEnhancedUpBBLinear::rec(double& fl, double& fr, size_t dim, grid_iterator& index) {
  size_t seq = index.seq();

  double fml = 0.0;
  double fmr = 0.0;
  fl = fr = 0.0;

  sgpp::base::level_t current_level;
  sgpp::base::index_t current_index;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec(fl, fml, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec(fmr, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }

  index.get(cur_algo_dim_, current_level, current_index);

  double fm = fml + fmr;
  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];

  // transposed operations:
  ptr_result_[(seq * this->numAlgoDims_) + dim] = fm;

  double tmp = (fm / 2.0) + (alpha_value / static_cast<double>(1 << (current_level + 1)));

  fl = tmp + fl;
  fr = tmp + fr;
}
#if 1
#if defined(__SSE3__)
void LaplaceEnhancedUpBBLinear::rec_LL(__m128d& fl, __m128d& fr, size_t dim, grid_iterator& index)
#else
void LaplaceEnhancedUpBBLinear::rec_LL(double& fl, double& fr, double& fl2, double& fr2,
                                       size_t dim, grid_iterator& index)
#endif
#else
#ifdef __SSE3__
void LaplaceEnhancedUpBBLinear::rec_LL(__m128d& fl, __m128d& fr, size_t dim, grid_iterator& index)
#else
void LaplaceEnhancedUpBBLinear::rec_LL(double& fl, double& fr, double& fl2, double& fr2,
                                       size_t dim, grid_iterator& index)
#endif
#endif
{
  size_t seq = index.seq();

#if 1
#if defined(__SSE3__)
  fl = _mm_set1_pd(0.0);
  fr = _mm_set1_pd(0.0);
  __m128d fml = _mm_set1_pd(0.0);
  __m128d fmr = _mm_set1_pd(0.0);
#else
  double fml = 0.0;
  double fmr = 0.0;
  double fml2 = 0.0;
  double fmr2 = 0.0;
  fl = fr = fl2 = fr2 = 0.0;
#endif
#else
#ifdef __SSE3__
  fl = _mm_set1_pd(0.0);
  fr = _mm_set1_pd(0.0);
  __m128d fml = _mm_set1_pd(0.0);
  __m128d fmr = _mm_set1_pd(0.0);
#else
  double fml = 0.0;
  double fmr = 0.0;
  double fml2 = 0.0;
  double fmr2 = 0.0;
  fl = fr = fl2 = fr2 = 0.0;
#endif
#endif
  sgpp::base::level_t current_level;
  sgpp::base::index_t current_index;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
#if 1
#if defined(__SSE3__)
      rec_LL(fl, fml, dim, index);
#else
      rec_LL(fl, fml, fl2, fml2, dim, index);
#endif
#else
#ifdef __SSE3__
      rec_LL(fl, fml, dim, index);
#else
      rec_LL(fl, fml, fl2, fml2, dim, index);
#endif
#endif
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
#if 1
#if defined(__SSE3__)
      rec_LL(fmr, fr, dim, index);
#else
      rec_LL(fmr, fr, fmr2, fr2, dim, index);
#endif
#else
#ifdef __SSE3__
      rec_LL(fmr, fr, dim, index);
#else
      rec_LL(fmr, fr, fmr2, fr2, dim, index);
#endif
#endif
    }

    index.up(cur_algo_dim_);
  }

  index.get(cur_algo_dim_, current_level, current_index);
#if 1
#if defined(__SSE3__)
  // mesh-width +1 level
  double h = 1.0 / static_cast<double>(1 << (current_level + 1));
  __m128d h_in = _mm_loaddup_pd(&h);
  __m128d fm = _mm_add_pd(fml, fmr);
  __m128d alpha = _mm_loadu_pd(&ptr_source_[(seq * this->numAlgoDims_) + dim]);
  _mm_storeu_pd(&ptr_result_[(seq * this->numAlgoDims_) + dim], fm);
  __m128d tmp = _mm_add_pd(_mm_mul_pd(fm, half_in_), _mm_mul_pd(alpha, h_in));
  fl = _mm_add_pd(fl, tmp);
  fr = _mm_add_pd(fr, tmp);
#else
  double fm = fml + fmr;
  double fm2 = fml2 + fmr2;
  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  double alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];

  // transposed operations:
  ptr_result_[(seq * this->numAlgoDims_) + dim] = fm;
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] = fm2;

  double tmp = (fm / 2.0) + (alpha_value / static_cast<double>(1 << (current_level + 1)));
  double tmp2 = (fm2 / 2.0) + (alpha_value2 / static_cast<double>(1 << (current_level + 1)));

  fl = tmp + fl;
  fr = tmp + fr;
  fl2 = tmp2 + fl2;
  fr2 = tmp2 + fr2;
#endif

#else
#ifdef __SSE3__
  // mesh-width +1 level
  double h = 1.0 / static_cast<double>(1 << (current_level + 1));
  __m128d h_in = _mm_loaddup_pd(&h);
  __m128d fm = _mm_add_pd(fml, fmr);
  __m128d alpha = _mm_loadu_pd(&ptr_source_[(seq * this->numAlgoDims_) + dim]);
  _mm_storeu_pd(&ptr_result_[(seq * this->numAlgoDims_) + dim], fm);
  __m128d tmp = _mm_add_pd(_mm_mul_pd(fm, half_in_), _mm_mul_pd(alpha, h_in));
  fl = _mm_add_pd(fl, tmp);
  fr = _mm_add_pd(fr, tmp);
#else
  double fm = fml + fmr;
  double fm2 = fml2 + fmr2;
  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  double alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];

  // transposed operations:
  ptr_result_[(seq * this->numAlgoDims_) + dim] = fm;
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] = fm2;

  double tmp = (fm / 2.0) + (alpha_value / static_cast<double>(1 << (current_level + 1)));
  double tmp2 = (fm2 / 2.0) + (alpha_value2 / static_cast<double>(1 << (current_level + 1)));

  fl = tmp + fl;
  fr = tmp + fr;
  fl2 = tmp2 + fl2;
  fr2 = tmp2 + fr2;
#endif
#endif
}

void LaplaceEnhancedUpBBLinear::rec_GL(double& fl, double& fr, size_t dim, grid_iterator& index) {
  size_t seq = index.seq();

  double fml = 0.0;
  double fmr = 0.0;
  fl = fr = 0.0;

  sgpp::base::level_t current_level;
  sgpp::base::index_t current_index;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec_GL(fl, fml, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec_GL(fmr, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }

  index.get(cur_algo_dim_, current_level, current_index);

  double fm = fml + fmr;
  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];

  // transposed operations:
  ptr_result_[(seq * this->numAlgoDims_) + dim] = 0.0;
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] = fm;

  double tmp = (fm / 2.0) + (alpha_value / static_cast<double>(1 << (current_level + 1)));

  fl = tmp + fl;
  fr = tmp + fr;
}

void LaplaceEnhancedUpBBLinear::rec_LG(double& fl, double& fr, size_t dim, grid_iterator& index) {
  size_t seq = index.seq();

  double fml = 0.0;
  double fmr = 0.0;
  fl = fr = 0.0;

  sgpp::base::level_t current_level;
  sgpp::base::index_t current_index;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec_LG(fl, fml, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      rec_LG(fmr, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }

  index.get(cur_algo_dim_, current_level, current_index);

  double fm = fml + fmr;
  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];

  // transposed operations:
  ptr_result_[(seq * this->numAlgoDims_) + dim] = fm;
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] = 0.0;

  double tmp = (fm / 2.0) + (alpha_value / static_cast<double>(1 << (current_level + 1)));

  fl = tmp + fl;
  fr = tmp + fr;
}

void LaplaceEnhancedUpBBLinear::rec_grad(size_t dim, grid_iterator& index) {
  size_t seq = index.seq();

  // Gradient
  ptr_result_[(seq * this->numAlgoDims_) + dim] = 0.0;

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

void LaplaceEnhancedUpBBLinear::recBB(double& fl, double& fr, size_t dim, grid_iterator& index) {
  size_t seq = index.seq();

  double fml = 0.0;
  double fmr = 0.0;
  fl = fr = 0.0;

  sgpp::base::level_t current_level;
  sgpp::base::index_t current_index;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB(fl, fml, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB(fmr, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }

  index.get(cur_algo_dim_, current_level, current_index);

  double fm = fml + fmr;
  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];

  // transposed operations:
  ptr_result_[(seq * this->numAlgoDims_) + dim] = fm;

  double tmp =
      ((fm / 2.0) + ((alpha_value / static_cast<double>(1 << (current_level + 1))) * q_));

  fl = tmp + fl;
  fr = tmp + fr;
}

void LaplaceEnhancedUpBBLinear::recBB_LL(double& fl, double& fr, double& fl2, double& fr2,
                                         size_t dim, grid_iterator& index) {
  size_t seq = index.seq();

  double fml = 0.0;
  double fmr = 0.0;
  double fml2 = 0.0;
  double fmr2 = 0.0;
  fl = fr = fl2 = fr2 = 0.0;

  sgpp::base::level_t current_level;
  sgpp::base::index_t current_index;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB_LL(fl, fml, fl2, fml2, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB_LL(fmr, fr, fmr2, fr2, dim, index);
    }

    index.up(cur_algo_dim_);
  }

  index.get(cur_algo_dim_, current_level, current_index);

  double fm = fml + fmr;
  double fm2 = fml2 + fmr2;
  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];
  double alpha_value2 = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];

  // transposed operations:
  ptr_result_[(seq * this->numAlgoDims_) + dim] = fm;
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] = fm2;

  double tmp = (fm / 2.0) + ((alpha_value / static_cast<double>(1 << (current_level + 1))) * q_);
  double tmp2 =
      (fm2 / 2.0) + ((alpha_value2 / static_cast<double>(1 << (current_level + 1))) * q_);

  fl = tmp + fl;
  fr = tmp + fr;
  fl2 = tmp2 + fl2;
  fr2 = tmp2 + fr2;
}

void LaplaceEnhancedUpBBLinear::recBB_GL(double& fl, double& fr, size_t dim,
                                         grid_iterator& index) {
  size_t seq = index.seq();

  double fml = 0.0;
  double fmr = 0.0;
  fl = fr = 0.0;

  sgpp::base::level_t current_level;
  sgpp::base::index_t current_index;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB_GL(fl, fml, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB_GL(fmr, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }

  index.get(cur_algo_dim_, current_level, current_index);

  double fm = fml + fmr;
  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim + 1];

  // transposed operations:
  ptr_result_[(seq * this->numAlgoDims_) + dim] = 0.0;
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] = fm;

  double tmp = (fm / 2.0) + ((alpha_value / static_cast<double>(1 << (current_level + 1))) * q_);

  fl = tmp + fl;
  fr = tmp + fr;
}

void LaplaceEnhancedUpBBLinear::recBB_LG(double& fl, double& fr, size_t dim,
                                         grid_iterator& index) {
  size_t seq = index.seq();

  double fml = 0.0;
  double fmr = 0.0;
  fl = fr = 0.0;

  sgpp::base::level_t current_level;
  sgpp::base::index_t current_index;

  if (!index.hint()) {
    index.leftChild(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB_LG(fl, fml, dim, index);
    }

    index.stepRight(cur_algo_dim_);

    if (!storage->isInvalidSequenceNumber(index.seq())) {
      recBB_LG(fmr, fr, dim, index);
    }

    index.up(cur_algo_dim_);
  }

  index.get(cur_algo_dim_, current_level, current_index);

  double fm = fml + fmr;
  double alpha_value = ptr_source_[(seq * this->numAlgoDims_) + dim];

  // transposed operations:
  ptr_result_[(seq * this->numAlgoDims_) + dim] = fm;
  ptr_result_[(seq * this->numAlgoDims_) + dim + 1] = 0.0;

  double tmp = (fm / 2.0) + ((alpha_value / static_cast<double>(1 << (current_level + 1))) * q_);

  fl = tmp + fl;
  fr = tmp + fr;
}

void LaplaceEnhancedUpBBLinear::recBB_grad(size_t dim, grid_iterator& index) {
  size_t seq = index.seq();

  // Gradient
  ptr_result_[(seq * this->numAlgoDims_) + dim] = 0.0;

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

}  // namespace pde
}  // namespace sgpp
