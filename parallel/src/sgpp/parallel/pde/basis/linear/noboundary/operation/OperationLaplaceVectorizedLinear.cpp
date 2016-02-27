// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <omp.h>

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <sgpp/parallel/pde/basis/linear/noboundary/operation/OperationLaplaceVectorizedLinear.hpp>

#include <sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/parallel/operation/HashGridStorageConverter.hpp>

#include <cmath>

#include <cstring>
#include <iostream>

#if defined(__SSE4_2__) || defined(__AVX__) || defined(__MIC__)
#include <immintrin.h>
#endif
#if defined(__FMA4__)
#include <x86intrin.h>
#endif

#ifdef __USEAVX128__
#undef __AVX__
#endif

#if defined(__MIC__)
#define VECTOR_SIZE 8
#define _mm512_and_pd(a, b) \
  _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#elif defined(__SSE4_2__) && !defined(__AVX__)
#define VECTOR_SIZE 2
#else
#define VECTOR_SIZE 4
#endif

#define REG_BCOUNT 2

#define BLOCK_LENGTH (REG_BCOUNT * VECTOR_SIZE)
#define PAGE_TRUE_CAPACITY 512

#define ROUND_DOWN(X, Y) (((X) / (Y)) * (Y))
#define ROUND_UP(X, Y) (std::ceil((1.0 * X) / (Y)) * (Y))

#include <sgpp/globaldef.hpp>
#include <algorithm>

namespace SGPP {
namespace parallel {

OperationLaplaceVectorizedLinear::OperationLaplaceVectorizedLinear(SGPP::base::GridStorage* storage)
    : storage(storage),
      level_(NULL),
      level_int_(NULL),
      index_(NULL),
      lcl_q_(NULL),
      lcl_q_inv_(NULL),
      lambda_(NULL),
      alpha_padded_(NULL),
      constants_(NULL),
      gradient_temp(NULL),
      l2dot_temp(NULL)
#if defined(STORE_MATRIX)
      ,
      operation_result_matrix_(NULL),
      operation_result_generated_(false)
#endif
{
  std::cout << "IN CONSTRUSTOR: OperationLaplaceVectorizedLinear" << std::endl;
  init_constants();
  init_grid_storage();

  this->lambda_ = new SGPP::base::DataVector(storage->getDimension());
  this->lambda_->setAll(1.0);
}

OperationLaplaceVectorizedLinear::OperationLaplaceVectorizedLinear(SGPP::base::GridStorage* storage,
                                                                   SGPP::base::DataVector& lambda)
    : storage(storage),
      level_(NULL),
      level_int_(NULL),
      index_(NULL),
      lcl_q_(NULL),
      lcl_q_inv_(NULL),
      lambda_(NULL),
      alpha_padded_(NULL),
      constants_(NULL),
      gradient_temp(NULL),
      l2dot_temp(NULL)
#if defined(STORE_MATRIX)
      ,
      operation_result_matrix_(NULL),
      operation_result_generated_(false)
#endif
{
  std::cout << "IN CONSTRUSTOR: OperationLaplaceVectorizedLinear" << std::endl;
  init_constants();
  init_grid_storage();

  this->lambda_ = new SGPP::base::DataVector(lambda);
}

OperationLaplaceVectorizedLinear::~OperationLaplaceVectorizedLinear() {
  double flop = static_cast<double>(26 * storage->getDimension() +
                                    storage->getDimension() * storage->getDimension()) *
                static_cast<double>(storage->getSize() * storage->getSize());
  double gflops = (all_iterations * flop / all_time) / 1000000000;
  double bandwidth = all_iterations *
                     static_cast<double>(sizeof(double) * storage->getSize() * storage->getSize()) /
                     all_time;
  std::cout << "IN OPERATOR : LAPLACE, GFLOPS :" << gflops
            << " BANDWIDTH :" << bandwidth / (1000000000.0) << " GB/s"
            << " ITERATIONS :" << all_iterations << " TIME :" << all_time << std::endl;
  // std::cout << gflops <<" " << bandwidth << " " << all_time << " " << all_iterations <<" "<<
  // storage->getSize() << '\n';
  delete this->level_;
  delete this->level_int_;
  delete this->index_;
  delete lcl_q_;
  delete lcl_q_inv_;

  if (this->lambda_ != NULL) delete this->lambda_;

  delete this->alpha_padded_;
  delete this->constants_;

#pragma omp parallel
  {
    delete gradient_temp[omp_get_thread_num()];
    delete l2dot_temp[omp_get_thread_num()];
  }
  delete gradient_temp;
  delete l2dot_temp;

#if defined(STORE_MATRIX)
  delete operation_result_matrix_;
#endif
}

void OperationLaplaceVectorizedLinear::reset() { init_grid_storage(); }

void OperationLaplaceVectorizedLinear::init_constants() {
  all_time = 0.0;
  all_iterations = 0.0;
#if defined(__SSE4_2__) || defined(__MIC__)
  this->constants_ = new SGPP::base::DataVector(0);

  this->constants_->append(0);
  this->constants_->append(0.5);
  this->constants_->append(2.0 / 3.0);
  this->constants_->append(1.0);
  this->constants_->append(2.0);

  uint64_t abs_mask = 0x7fffffffffffffff;
  this->constants_->append(*(reinterpret_cast<double*>(&abs_mask)));
#endif

  process_index = 0;
  process_count = 1;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &process_index);
  MPI_Comm_size(MPI_COMM_WORLD, &process_count);
#endif
}

void OperationLaplaceVectorizedLinear::init_grid_storage() {
  if (this->level_) delete this->level_;

  this->level_ = new SGPP::base::DataMatrix(storage->getSize(), storage->getDimension());

  if (this->level_int_) delete this->level_int_;

  this->level_int_ = new SGPP::base::DataMatrix(storage->getSize(), storage->getDimension());

  if (this->index_) delete this->index_;

  this->index_ = new SGPP::base::DataMatrix(storage->getSize(), storage->getDimension());

  if (this->lcl_q_) delete this->lcl_q_;

  lcl_q_ = new SGPP::base::DataVector(this->storage->getDimension());
  double* lcl_q_ptr_ = lcl_q_->getPointer();

  if (this->lcl_q_inv_) delete this->lcl_q_inv_;

  lcl_q_inv_ = new SGPP::base::DataVector(this->storage->getDimension());
  double* lcl_q_inv_ptr_ = lcl_q_inv_->getPointer();

#if defined(__MIC__)
  SGPP::parallel::HashGridStorageConverter::getLevelIndexArraysForEvalTLBOptimized(
      storage, *(this->level_), *(this->index_), SGPP::parallel::MIC, BLOCK_LENGTH);
  SGPP::parallel::HashGridStorageConverter::getLevelForIntegralTLBOptimized(
      storage, *(this->level_int_), SGPP::parallel::MIC, BLOCK_LENGTH);
#elif defined(__SSE4_2__) || defined(__AVX__)
  SGPP::parallel::HashGridStorageConverter::getLevelIndexArraysForEvalTLBOptimized(
      storage, *(this->level_), *(this->index_), SGPP::parallel::X86SIMD, BLOCK_LENGTH);
  SGPP::parallel::HashGridStorageConverter::getLevelForIntegralTLBOptimized(
      storage, *(this->level_int_), SGPP::parallel::X86SIMD, BLOCK_LENGTH);
#else
  storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
  storage->getLevelForIntegral(*(this->level_int_));
#endif

#if defined(__MIC__) || defined(__SSE4_2__) || defined(__AVX__)
  std::size_t padded_size = this->level_->getNcols();
#else
  std::size_t padded_size = this->level_->getNrows();
#endif

  if (this->alpha_padded_) delete this->alpha_padded_;

  this->alpha_padded_ = new SGPP::base::DataVector(padded_size);
  this->alpha_padded_->setAll(0.0);

  size_t single_process_portion = (this->storage->getSize() / process_count) + 1;

  all_i_start.clear();
  all_i_size.clear();
  send_start.clear();
  send_size.clear();
  recv_start.clear();
  recv_size.clear();

  for (int i = 0; i < process_count; ++i) {
    int process_start = static_cast<int>(i * single_process_portion);
    all_i_start.push_back(process_start);
    int process_portion =
        (i == process_count - 1)
            ? static_cast<int>(this->storage->getSize() - i * single_process_portion)
            : static_cast<int>(single_process_portion);

    process_portion = std::min<int>(
        process_portion, static_cast<int>(this->storage->getSize() - i * single_process_portion));
    process_portion = std::max<int>(process_portion, 0);

    all_i_size.push_back(process_portion);
  }

  for (int i = 0; i < process_count; ++i) {
    send_start.push_back(static_cast<int>(single_process_portion * process_index));
    int process_send_size =
        (process_index == process_count - 1)
            ? static_cast<int>(this->storage->getSize() - single_process_portion * process_index)
            : static_cast<int>(single_process_portion);

    process_send_size = std::min<int>(
        process_send_size,
        static_cast<int>(this->storage->getSize() - single_process_portion * process_index));
    process_send_size = std::max<int>(process_send_size, 0);
    // process_send_size = (i == process_index)? 0: process_send_size;

    send_size.push_back(process_send_size);

    recv_start.push_back(static_cast<int>(single_process_portion * i));
    int process_recv_size =
        (i == process_count - 1)
            ? static_cast<int>(this->storage->getSize() - single_process_portion * i)
            : static_cast<int>(single_process_portion);

    process_recv_size = std::min<int>(
        process_recv_size, static_cast<int>(this->storage->getSize() - single_process_portion * i));
    process_recv_size = std::max<int>(process_recv_size, 0);
    // process_recv_size = (i == process_index)? 0 : process_recv_size;

    recv_size.push_back(process_recv_size);
  }

#pragma omp parallel
  {
    if (this->gradient_temp) delete this->gradient_temp[omp_get_thread_num()];

    if (this->l2dot_temp) delete this->l2dot_temp[omp_get_thread_num()];

#pragma omp barrier

#pragma omp single
    {
      if (this->gradient_temp) delete this->gradient_temp;

      if (this->l2dot_temp) delete this->l2dot_temp;

      this->gradient_temp = new SGPP::base::DataVector*[omp_get_num_threads()];
      this->l2dot_temp = new SGPP::base::DataVector*[omp_get_num_threads()];
    }
#pragma omp barrier

    // std::cout << "OMP THREAD :" << omp_get_thread_num() << std::endl;

    gradient_temp[omp_get_thread_num()] =
        new SGPP::base::DataVector(VECTOR_SIZE * this->storage->getDimension() * REG_BCOUNT);
    l2dot_temp[omp_get_thread_num()] =
        new SGPP::base::DataVector(VECTOR_SIZE * this->storage->getDimension() * REG_BCOUNT);
  }

  // fill q array
  for (size_t d = 0; d < this->storage->getDimension(); d++) {
    SGPP::base::BoundingBox* boundingBox = this->storage->getBoundingBox();
    lcl_q_ptr_[d] = boundingBox->getIntervalWidth(d);
    lcl_q_inv_ptr_[d] = 1.0 / boundingBox->getIntervalWidth(d);
  }

#if defined(STORE_MATRIX)
  size_t result_matrix_rows = all_i_size[process_index];
#if defined(__MIC__) || defined(__SSE4_2__) || defined(__AVX__)
  size_t result_matrix_cols = this->level_->getNcols();
#else
  size_t result_matrix_cols = this->level_->getNrows();
#endif
  // check if matrix fits in memory
  char* matrix_max_size_gb_str = getenv("SGPP_PDE_MATRIX_SIZE_GB");
  size_t matrix_max_size_gb = 0;

  if (matrix_max_size_gb_str == NULL || strlen(matrix_max_size_gb_str) == 0) {
    matrix_max_size_gb = 4;
  } else {
    std::stringstream temp_stream(matrix_max_size_gb_str);
    temp_stream >> matrix_max_size_gb;
  }

  size_t matrix_needed_size_bytes = result_matrix_rows * result_matrix_cols * sizeof(double);

  if (matrix_needed_size_bytes > matrix_max_size_gb * 1024 * 1024 * 1024) {
    size_t matrix_needed_size_gb = static_cast<size_t>(
        ceil(static_cast<double>(matrix_needed_size_bytes) / (1024 * 1024 * 1024)));
    char exception_string[512];
    snprintf(exception_string, sizeof(exception_string)
            "OperationLaplaceVectorizedLinear::init : More memory (= %i GB) needed to store the "
            "operation matrix, Please set the SGPP_PDE_MATRIX_SIZE_GB environment variable!",
            static_cast<int>(matrix_needed_size_gb));

    std::cerr << exception_string << std::endl;
    throw SGPP::base::operation_exception(exception_string);
  }

  if (operation_result_matrix_) delete operation_result_matrix_;

  operation_result_matrix_ = new SGPP::base::DataMatrix(result_matrix_rows, result_matrix_cols);
  operation_result_generated_ = false;

#pragma omp parallel
  {
    size_t padded_size = this->operation_result_matrix_->getNcols();
    size_t thr_start;
    size_t thr_end;
    SGPP::parallel::PartitioningTool::getOpenMPPartitionSegment(0, result_matrix_rows, &thr_start,
                                                                &thr_end);

    for (size_t i = thr_start; i < thr_end; i++) {
      double* operation_result_dest_ptr =
          operation_result_matrix_->getPointer() + (i)*operation_result_matrix_->getNcols();

      for (size_t j = 0; j < padded_size; ++j) {
        operation_result_dest_ptr[j] = 0.0;
      }
    }
  }
#endif
}

double OperationLaplaceVectorizedLinear::gradient(size_t i, size_t j, size_t dim) {
  double grad;

  double i_level_grad = level_->get(i, dim);
  double i_index_grad = index_->get(i, dim);
  double j_level_grad = level_->get(j, dim);
  double j_index_grad = index_->get(j, dim);

  // only affects the diagonal of the stiffness matrix

  bool doGrad = ((i_level_grad == j_level_grad) && (i_index_grad == j_index_grad));
  grad = i_level_grad * 2.0 * static_cast<double>(static_cast<int>(doGrad));

  return (grad * *(lcl_q_inv_->getPointer() + dim));
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
  // double l1d = lid*(selector) + ljd*(!selector);
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
  double res_two = temp_res * overlap;  // Now mask result

  return (res_one * (lid == ljd) + res_two * (lid != ljd)) * *(lcl_q_->getPointer() + dim);
}

#if defined(__MIC__)
void mic_mult(size_t process_i_start, size_t process_i_end, SGPP::base::DataVector& result) {
  {
    std::size_t padded_size = this->level_->getNcols();
    double* constants = this->constants_->getPointer();  // {0, 0.5, 2.0 / 3.0, 1, 2};
    double* result_ptr_ = result.getPointer();
    double* level_ptr_ = this->level_->getPointer();
    double* level_int_ptr_ = this->level_int_->getPointer();
    double* index_ptr_ = this->index_->getPointer();
    double* gradient_temp_ptr = gradient_temp[omp_get_thread_num()]->getPointer();
    double* l2dot_temp_ptr = l2dot_temp[omp_get_thread_num()]->getPointer();
    double* lambda_ptr_ = this->lambda_->getPointer();
    size_t temp_cols = this->storage->getDimension() * VECTOR_SIZE;
    double* lcl_q_temp_ptr_ = lcl_q_->getPointer();
    double* lcl_q_inv_temp_ptr_ = lcl_q_inv_->getPointer();
    double* alpha_padded_temp_ptr_ = alpha_padded_->getPointer();
    size_t max_dims = this->storage->getDimension();
    size_t page_cap_rounded = max_dims * BLOCK_LENGTH;
    __m512d mm_zero =
        _mm512_extload_pd(constants + 0, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
    __m512d mm_half =
        _mm512_extload_pd(constants + 1, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
    __m512d mm_two_thirds =
        _mm512_extload_pd(constants + 2, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
    __m512d mm_one =
        _mm512_extload_pd(constants + 3, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
    __m512d mm_two =
        _mm512_extload_pd(constants + 4, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
    __m512d mm_abs =
        _mm512_extload_pd(constants + 5, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
    size_t thr_start;
    size_t thr_end;
    SGPP::parallel::PartitioningTool::getOpenMPPartitionSegment(process_i_start, process_i_end,
                                                                &thr_start, &thr_end);
    for (size_t i = thr_start; i < thr_end; i++) {
      __m512d mm_result = mm_zero;
      __m512d mm_result2 = mm_zero;
      size_t i_page = i / BLOCK_LENGTH;
      size_t i_offset = i_page * page_cap_rounded + i % BLOCK_LENGTH;
      double* temp_level_ptr = level_ptr_;
      double* temp_level_int_ptr = level_int_ptr_;
      double* temp_index_ptr = index_ptr_;
      for (size_t j = 0; j < padded_size; j += VECTOR_SIZE * REG_BCOUNT) {
        double* gradient_temp_ptr1 = gradient_temp_ptr;
        double* l2dot_temp_ptr1 = l2dot_temp_ptr;
        size_t i_idx = i_offset;
        for (size_t dim = 0; dim < max_dims; dim++) {
          __m512d mm_lcl_q = _mm512_extload_pd(lcl_q_temp_ptr_ + dim, _MM_UPCONV_PD_NONE,
                                               _MM_BROADCAST_1X8, _MM_HINT_NONE);
          __m512d mm_lcl_q_inv = _mm512_extload_pd(lcl_q_inv_temp_ptr_ + dim, _MM_UPCONV_PD_NONE,
                                                   _MM_BROADCAST_1X8, _MM_HINT_NONE);
          __m512d mm_lid = _mm512_extload_pd(level_ptr_ + i_idx, _MM_UPCONV_PD_NONE,
                                             _MM_BROADCAST_1X8, _MM_HINT_NONE);
          __m512d mm_iid = _mm512_extload_pd(index_ptr_ + i_idx, _MM_UPCONV_PD_NONE,
                                             _MM_BROADCAST_1X8, _MM_HINT_NONE);
          __m512d mm_ljd = _mm512_load_pd(temp_level_ptr);
          __m512d mm_ijd = _mm512_load_pd(temp_index_ptr);
          __mmask8 mm_doGrad =
              (__mmask8)(_mm512_kand(_mm512_cmp_pd_mask(mm_lid, mm_ljd, _MM_CMPINT_EQ),
                                     _mm512_cmp_pd_mask(mm_iid, mm_ijd, _MM_CMPINT_EQ)));  // +3
          __m512d mm_grad = _mm512_mask_mul_pd(mm_zero, mm_doGrad, mm_lid, mm_two);
          mm_grad = _mm512_mul_pd(mm_grad, mm_lcl_q_inv);  // 1
          _mm512_store_pd(gradient_temp_ptr1, mm_grad);
          __m512d mm_in_lid = _mm512_extload_pd(level_int_ptr_ + i_idx, _MM_UPCONV_PD_NONE,
                                                _MM_BROADCAST_1X8, _MM_HINT_NONE);
          __m512d mm_in_ljd = _mm512_load_pd(temp_level_int_ptr);
          __m512d mm_res_one =
              _mm512_mask_mul_pd(mm_zero, _mm512_cmp_pd_mask(mm_iid, mm_ijd, _MM_CMPINT_EQ),
                                 mm_two_thirds, mm_in_lid);                          // 1+2
          __mmask8 mm_selector = _mm512_cmp_pd_mask(mm_lid, mm_ljd, _MM_CMPINT_LE);  // +6
          __m512d mm_i1d = _mm512_mask_blend_pd(mm_selector, mm_iid, mm_ijd);
          __m512d mm_in_l1d = _mm512_mask_blend_pd(mm_selector, mm_in_lid, mm_in_ljd);
          __m512d mm_in_l2d = _mm512_mask_blend_pd(mm_selector, mm_in_ljd, mm_in_lid);
          __m512d mm_i2d = _mm512_mask_blend_pd(mm_selector, mm_ijd, mm_iid);
          __m512d mm_l2d = _mm512_mask_blend_pd(mm_selector, mm_ljd, mm_lid);
          __m512d mm_q = _mm512_mul_pd(_mm512_sub_pd(mm_i1d, mm_one), mm_in_l1d);  // 2 flop
          __m512d mm_p = _mm512_mul_pd(_mm512_add_pd(mm_i1d, mm_one), mm_in_l1d);  // 2 flop
          __mmask8 mm_overlap = _mm512_cmp_pd_mask(
              _mm512_max_pd(mm_q, _mm512_mul_pd(_mm512_sub_pd(mm_i2d, mm_one), mm_in_l2d)),
              _mm512_min_pd(mm_p, _mm512_mul_pd(_mm512_add_pd(mm_i2d, mm_one), mm_in_l2d)),
              _MM_CMPINT_LT);  // 6+1
          __m512d mm_temp_res = _mm512_sub_pd(
              _mm512_sub_pd(
                  mm_two,
                  _mm512_and_pd(mm_abs, (_mm512_sub_pd(_mm512_mul_pd(mm_l2d, mm_q), mm_i2d)))),
              _mm512_and_pd(mm_abs, (_mm512_sub_pd(_mm512_mul_pd(mm_l2d, mm_p),
                                                   mm_i2d))));  // 8 flops
          __m512d mm_res_two = _mm512_mask_mul_pd(mm_zero, mm_overlap, mm_temp_res,
                                                  _mm512_mul_pd(mm_half, mm_in_l1d));
          mm_selector = _mm512_cmp_pd_mask(mm_lid, mm_ljd, _MM_CMPINT_NE);  // +1
          __m512d mm_val = _mm512_mask_blend_pd(mm_selector, mm_res_one, mm_res_two);
          mm_val = _mm512_mul_pd(mm_val, mm_lcl_q);  // 1 flop
          _mm512_store_pd(l2dot_temp_ptr1, mm_val);
          ////////////////////////////////////////////////////////
          __m512d mm_ljd2 = _mm512_load_pd(temp_level_ptr + VECTOR_SIZE);
          __m512d mm_ijd2 = _mm512_load_pd(temp_index_ptr + VECTOR_SIZE);
          __mmask8 mm_doGrad2 = (__mmask8)(
              _mm512_kand(_mm512_cmp_pd_mask(mm_lid, mm_ljd2, _MM_CMPINT_EQ),
                          _mm512_cmp_pd_mask(mm_iid, mm_ijd2, _MM_CMPINT_EQ)));  // 1 // +2
          __m512d mm_grad2 = _mm512_mask_mul_pd(mm_zero, mm_doGrad2, mm_lid, mm_two);
          mm_grad2 = _mm512_mul_pd(mm_grad2, mm_lcl_q_inv);  // 1
          _mm512_store_pd(gradient_temp_ptr1 + temp_cols, mm_grad2);
          __m512d mm_in_ljd2 = _mm512_load_pd(temp_level_int_ptr + VECTOR_SIZE);
          __m512d mm_res_one2 =
              _mm512_mask_mul_pd(mm_zero, _mm512_cmp_pd_mask(mm_iid, mm_ijd2, _MM_CMPINT_EQ),
                                 mm_two_thirds, mm_in_lid);  // 1+2
          __mmask8 mm_selector2 = _mm512_cmp_pd_mask(mm_lid, mm_ljd2, _MM_CMPINT_LE);
          __m512d mm_i1d2 = _mm512_mask_blend_pd(mm_selector2, mm_iid, mm_ijd2);
          __m512d mm_in_l1d2 = _mm512_mask_blend_pd(mm_selector2, mm_in_lid, mm_in_ljd2);
          __m512d mm_in_l2d2 = _mm512_mask_blend_pd(mm_selector2, mm_in_ljd2, mm_in_lid);
          __m512d mm_i2d2 = _mm512_mask_blend_pd(mm_selector2, mm_ijd2, mm_iid);
          __m512d mm_l2d2 = _mm512_mask_blend_pd(mm_selector2, mm_ljd2, mm_lid);
          __m512d mm_q2 = _mm512_mul_pd(_mm512_sub_pd(mm_i1d2, mm_one), mm_in_l1d2);  // 2 flop
          __m512d mm_p2 = _mm512_mul_pd(_mm512_add_pd(mm_i1d2, mm_one), mm_in_l1d2);  // 2 flop
          __mmask8 mm_overlap2 = _mm512_cmp_pd_mask(
              _mm512_max_pd(mm_q2, _mm512_mul_pd(_mm512_sub_pd(mm_i2d2, mm_one), mm_in_l2d2)),
              _mm512_min_pd(mm_p2, _mm512_mul_pd(_mm512_add_pd(mm_i2d2, mm_one), mm_in_l2d2)),
              _MM_CMPINT_LT);  // 6 flop // +1
          __m512d mm_temp_res2 = _mm512_sub_pd(
              _mm512_sub_pd(
                  mm_two,
                  _mm512_and_pd(mm_abs, (_mm512_sub_pd(_mm512_mul_pd(mm_l2d2, mm_q2), mm_i2d2)))),
              _mm512_and_pd(mm_abs, (_mm512_sub_pd(_mm512_mul_pd(mm_l2d2, mm_p2),
                                                   mm_i2d2))));  // 8 flops
          __m512d mm_res_two2 = _mm512_mask_mul_pd(mm_zero, mm_overlap2, mm_temp_res2,
                                                   _mm512_mul_pd(mm_half, mm_in_l1d2));
          mm_selector2 = _mm512_cmp_pd_mask(mm_lid, mm_ljd2, _MM_CMPINT_NE);
          __m512d mm_val2 = _mm512_mask_blend_pd(mm_selector2, mm_res_one2, mm_res_two2);
          mm_val2 = _mm512_mul_pd(mm_val2, mm_lcl_q);  // 1 flop
          _mm512_store_pd(l2dot_temp_ptr1 + temp_cols, mm_val2);
          gradient_temp_ptr1 += VECTOR_SIZE;
          l2dot_temp_ptr1 += VECTOR_SIZE;
          temp_level_ptr += BLOCK_LENGTH;
          temp_level_int_ptr += BLOCK_LENGTH;
          temp_index_ptr += BLOCK_LENGTH;
          i_idx += BLOCK_LENGTH;
        }
        gradient_temp_ptr1 = gradient_temp_ptr;
        for (size_t d_outer = 0; d_outer < max_dims; d_outer++) {
          // D + 2
          __m512d mm_element = _mm512_load_pd(alpha_padded_temp_ptr_ + j);
          __m512d mm_temp = _mm512_load_pd(gradient_temp_ptr1);
          __m512d mm_element2 = _mm512_load_pd(alpha_padded_temp_ptr_ + j + VECTOR_SIZE);
          __m512d mm_temp2 = _mm512_load_pd(gradient_temp_ptr1 + temp_cols);
          mm_element = _mm512_mul_pd(mm_element, mm_temp);
          mm_element2 = _mm512_mul_pd(mm_element2, mm_temp2);
          l2dot_temp_ptr1 = l2dot_temp_ptr;
          for (size_t d_inner = 0; d_inner < d_outer; d_inner++) {
            // element *= ((l2dot(i, j, d_inner)*(d_outer != d_inner)) + (gradient(i, j,
            // d_inner)*(d_outer == d_inner)));
            __m512d mm_temp = _mm512_load_pd(l2dot_temp_ptr1);
            __m512d mm_temp2 = _mm512_load_pd(l2dot_temp_ptr1 + temp_cols);
            mm_element = _mm512_mul_pd(mm_element, mm_temp);
            mm_element2 = _mm512_mul_pd(mm_element2, mm_temp2);
            l2dot_temp_ptr1 += VECTOR_SIZE;
          }
          l2dot_temp_ptr1 += VECTOR_SIZE;
          for (size_t d_inner = d_outer + 1; d_inner < max_dims; d_inner++) {
            // element *= ((l2dot(i, j, d_inner)*(d_outer != d_inner)) + (gradient(i, j,
            // d_inner)*(d_outer == d_inner)));
            __m512d mm_temp = _mm512_load_pd(l2dot_temp_ptr1);
            __m512d mm_temp2 = _mm512_load_pd(l2dot_temp_ptr1 + temp_cols);
            mm_element = _mm512_mul_pd(mm_element, mm_temp);
            mm_element2 = _mm512_mul_pd(mm_element2, mm_temp2);
            l2dot_temp_ptr1 += VECTOR_SIZE;
          }
          __m512d mm_lambda = _mm512_extload_pd(lambda_ptr_ + d_outer, _MM_UPCONV_PD_NONE,
                                                _MM_BROADCAST_1X8, _MM_HINT_NONE);
          mm_element = _mm512_mul_pd(mm_element, mm_lambda);
          mm_element2 = _mm512_mul_pd(mm_element2, mm_lambda);
          mm_result = _mm512_add_pd(mm_result, mm_element);
          mm_result2 = _mm512_add_pd(mm_result2, mm_element2);
          gradient_temp_ptr1 += VECTOR_SIZE;
        }
      }
      mm_result = _mm512_add_pd(mm_result, mm_result2);
      result_ptr_[i] += _mm512_reduce_add_pd(mm_result);
    }
  }
}
#endif

void OperationLaplaceVectorizedLinear::mult(SGPP::base::DataVector& alpha,
                                            SGPP::base::DataVector& result) {
  result.setAll(0.0);

  stopWatch.start();
  size_t process_i_start = all_i_start[process_index];
  size_t process_i_end = process_i_start + all_i_size[process_index];

// std::cout << "PROCESS :" << process_index << " START :" << process_i_start << " END :" <<
// process_i_end << " COUNT :" << (process_i_end - process_i_start) <<std::endl;

#if defined(STORE_MATRIX)

  if (!operation_result_generated_) {
    operation_result_generated_ = true;

    alpha_padded_->setAll(1.0);
#else
  std::size_t original_size = alpha.getSize();
  memcpy(alpha_padded_->getPointer(), alpha.getPointer(), original_size * sizeof(double));
#endif

#if defined(__MIC__)

#pragma omp parallel
    mic_mult(process_i_start, process_i_end, result);

#elif defined(__SSE4_2__) && defined(__AVX__)

#pragma omp parallel
  {
    std::size_t padded_size = this->level_->getNcols();
    double* constants = this->constants_->getPointer();  // {0, 0.5, 2.0 / 3.0, 1, 2};
#if !defined(STORE_MATRIX)
    double* result_ptr_ = result.getPointer();
#endif
    double* level_ptr_ = this->level_->getPointer();
    double* level_int_ptr_ = this->level_int_->getPointer();
    double* index_ptr_ = this->index_->getPointer();

    double* gradient_temp_ptr = gradient_temp[omp_get_thread_num()]->getPointer();
    double* l2dot_temp_ptr = l2dot_temp[omp_get_thread_num()]->getPointer();
    double* lambda_ptr_ = this->lambda_->getPointer();
    size_t temp_cols = this->storage->getDimension() * VECTOR_SIZE;

    double* lcl_q_temp_ptr_ = lcl_q_->getPointer();
    double* lcl_q_inv_temp_ptr_ = lcl_q_inv_->getPointer();
    double* alpha_padded_temp_ptr_ = alpha_padded_->getPointer();

    size_t max_dims = this->storage->getDimension();
    size_t page_cap_rounded = max_dims * BLOCK_LENGTH;

    __m256d mm_half = _mm256_broadcast_sd(constants + 1);
    __m256d mm_two_thirds = _mm256_broadcast_sd(constants + 2);
    __m256d mm_one = _mm256_broadcast_sd(constants + 3);
    __m256d mm_two = _mm256_broadcast_sd(constants + 4);
    __m256d mm_abs = _mm256_broadcast_sd(constants + 5);

    size_t thr_start;
    size_t thr_end;
    SGPP::parallel::PartitioningTool::getOpenMPPartitionSegment(process_i_start, process_i_end,
                                                                &thr_start, &thr_end);

    for (size_t i = thr_start; i < thr_end; i++) {
      __m256d mm_result = _mm256_setzero_pd();
      __m256d mm_result2 = _mm256_setzero_pd();

      size_t i_page = i / BLOCK_LENGTH;
      size_t i_offset = i_page * page_cap_rounded + i % BLOCK_LENGTH;

      double* temp_level_ptr = level_ptr_;
      double* temp_level_int_ptr = level_int_ptr_;
      double* temp_index_ptr = index_ptr_;

      for (size_t j = 0; j < padded_size; j += VECTOR_SIZE * REG_BCOUNT) {
#if defined(STORE_MATRIX)
        mm_result = _mm256_setzero_pd();
        mm_result2 = _mm256_setzero_pd();
#endif
        double* gradient_temp_ptr1 = gradient_temp_ptr;
        double* l2dot_temp_ptr1 = l2dot_temp_ptr;

        size_t i_idx = i_offset;

        for (size_t dim = 0; dim < max_dims; dim++) {
          __m256d mm_lcl_q = _mm256_broadcast_sd(lcl_q_temp_ptr_ + dim);
          __m256d mm_lcl_q_inv = _mm256_broadcast_sd(lcl_q_inv_temp_ptr_ + dim);

          __m256d mm_lid = _mm256_broadcast_sd(level_ptr_ + i_idx);
          __m256d mm_iid = _mm256_broadcast_sd(index_ptr_ + i_idx);
          __m256d mm_ljd = _mm256_load_pd(temp_level_ptr);
          __m256d mm_ijd = _mm256_load_pd(temp_index_ptr);

          __m256d mm_doGrad = _mm256_and_pd(_mm256_cmp_pd(mm_lid, mm_ljd, _CMP_EQ_OQ),
                                            _mm256_cmp_pd(mm_iid, mm_ijd, _CMP_EQ_OQ));  // +3

          __m256d mm_grad = _mm256_mul_pd(mm_lid, _mm256_and_pd(mm_two, mm_doGrad));  // 1+1

          mm_grad = _mm256_mul_pd(mm_grad, mm_lcl_q_inv);  // 1

          _mm256_store_pd(gradient_temp_ptr1, mm_grad);

          __m256d mm_in_lid = _mm256_broadcast_sd(level_int_ptr_ + i_idx);
          __m256d mm_in_ljd = _mm256_load_pd(temp_level_int_ptr);

          __m256d mm_res_one = _mm256_mul_pd(
              mm_two_thirds,
              _mm256_and_pd(mm_in_lid, _mm256_cmp_pd(mm_iid, mm_ijd, _CMP_EQ_OQ)));  // 1+2

          __m256d mm_selector = _mm256_cmp_pd(mm_lid, mm_ljd, _CMP_LE_OQ);  // +6
          __m256d mm_i1d = _mm256_blendv_pd(mm_iid, mm_ijd, mm_selector);
          __m256d mm_in_l1d = _mm256_blendv_pd(mm_in_lid, mm_in_ljd, mm_selector);
          __m256d mm_in_l2d = _mm256_blendv_pd(mm_in_ljd, mm_in_lid, mm_selector);
          __m256d mm_i2d = _mm256_blendv_pd(mm_ijd, mm_iid, mm_selector);
          __m256d mm_l2d = _mm256_blendv_pd(mm_ljd, mm_lid, mm_selector);

          __m256d mm_q = _mm256_mul_pd(_mm256_sub_pd(mm_i1d, mm_one), mm_in_l1d);  // 2 flop
          __m256d mm_p = _mm256_mul_pd(_mm256_add_pd(mm_i1d, mm_one), mm_in_l1d);  // 2 flop
          __m256d mm_overlap = _mm256_cmp_pd(
              _mm256_max_pd(mm_q, _mm256_mul_pd(_mm256_sub_pd(mm_i2d, mm_one), mm_in_l2d)),
              _mm256_min_pd(mm_p, _mm256_mul_pd(_mm256_add_pd(mm_i2d, mm_one), mm_in_l2d)),
              _CMP_LT_OQ);  // 6+1

          __m256d mm_temp_res = _mm256_sub_pd(
              _mm256_sub_pd(
                  mm_two,
                  _mm256_and_pd(mm_abs, (_mm256_sub_pd(_mm256_mul_pd(mm_l2d, mm_q), mm_i2d)))),
              _mm256_and_pd(mm_abs, (_mm256_sub_pd(_mm256_mul_pd(mm_l2d, mm_p),
                                                   mm_i2d))));                          // 8 flops
          mm_temp_res = _mm256_mul_pd(mm_temp_res, _mm256_mul_pd(mm_half, mm_in_l1d));  // 2 flops
          __m256d mm_res_two = _mm256_and_pd(mm_temp_res, mm_overlap);  // Now mask result //+1
          mm_selector = _mm256_cmp_pd(mm_lid, mm_ljd, _CMP_NEQ_OQ);     // +1

          __m256d mm_val = _mm256_blendv_pd(mm_res_one, mm_res_two, mm_selector);  // +1
          mm_val = _mm256_mul_pd(mm_val, mm_lcl_q);                                // 1 flop
          _mm256_store_pd(l2dot_temp_ptr1, mm_val);

          ////////////////////////////////////////////////////////
          __m256d mm_ljd2 = _mm256_load_pd(temp_level_ptr + VECTOR_SIZE);
          __m256d mm_ijd2 = _mm256_load_pd(temp_index_ptr + VECTOR_SIZE);
          __m256d mm_doGrad2 =
              _mm256_and_pd(_mm256_cmp_pd(mm_lid, mm_ljd2, _CMP_EQ_OQ),
                            _mm256_cmp_pd(mm_iid, mm_ijd2, _CMP_EQ_OQ));  // 1 // +2

          __m256d mm_grad2 = _mm256_mul_pd(mm_lid, _mm256_and_pd(mm_two, mm_doGrad2));  // 2

          mm_grad2 = _mm256_mul_pd(mm_grad2, mm_lcl_q_inv);  // 1
          _mm256_store_pd(gradient_temp_ptr1 + temp_cols, mm_grad2);

          __m256d mm_in_ljd2 = _mm256_load_pd(temp_level_int_ptr + VECTOR_SIZE);

          __m256d mm_res_one2 = _mm256_mul_pd(
              mm_two_thirds,
              _mm256_and_pd(mm_in_lid, _mm256_cmp_pd(mm_iid, mm_ijd2, _CMP_EQ_OQ)));  // 2 // +1
          __m256d mm_selector2 = _mm256_cmp_pd(mm_lid, mm_ljd2, _CMP_LE_OQ);
          __m256d mm_i1d2 = _mm256_blendv_pd(mm_iid, mm_ijd2, mm_selector2);
          __m256d mm_in_l1d2 = _mm256_blendv_pd(mm_in_lid, mm_in_ljd2, mm_selector2);
          __m256d mm_i2d2 = _mm256_blendv_pd(mm_ijd2, mm_iid, mm_selector2);
          __m256d mm_l2d2 = _mm256_blendv_pd(mm_ljd2, mm_lid, mm_selector2);
          __m256d mm_in_l2d2 = _mm256_blendv_pd(mm_in_ljd2, mm_in_lid, mm_selector2);

          __m256d mm_q2 = _mm256_mul_pd(_mm256_sub_pd(mm_i1d2, mm_one), mm_in_l1d2);  // 2 flop
          __m256d mm_p2 = _mm256_mul_pd(_mm256_add_pd(mm_i1d2, mm_one), mm_in_l1d2);  // 2 flop

          __m256d mm_overlap2 = _mm256_cmp_pd(
              _mm256_max_pd(mm_q2, _mm256_mul_pd(_mm256_sub_pd(mm_i2d2, mm_one), mm_in_l2d2)),
              _mm256_min_pd(mm_p2, _mm256_mul_pd(_mm256_add_pd(mm_i2d2, mm_one), mm_in_l2d2)),
              _CMP_LT_OQ);  // 6 flop // +1

          __m256d mm_temp_res2 = _mm256_sub_pd(
              _mm256_sub_pd(
                  mm_two,
                  _mm256_and_pd(mm_abs, (_mm256_sub_pd(_mm256_mul_pd(mm_l2d2, mm_q2), mm_i2d2)))),
              _mm256_and_pd(mm_abs, (_mm256_sub_pd(_mm256_mul_pd(mm_l2d2, mm_p2),
                                                   mm_i2d2))));  // 8 flops

          mm_temp_res2 = _mm256_mul_pd(mm_temp_res2, _mm256_mul_pd(mm_half,
                                                                   mm_in_l1d2));  // 2 flops
          __m256d mm_res_two2 = _mm256_and_pd(mm_temp_res2,
                                              mm_overlap2);  // Now mask result //1 flop

          mm_selector2 = _mm256_cmp_pd(mm_lid, mm_ljd2, _CMP_NEQ_OQ);

          __m256d mm_val2 = _mm256_blendv_pd(mm_res_one2, mm_res_two2, mm_selector2);

          mm_val2 = _mm256_mul_pd(mm_val2, mm_lcl_q);  // 1 flop

          _mm256_store_pd(l2dot_temp_ptr1 + temp_cols, mm_val2);

          gradient_temp_ptr1 += VECTOR_SIZE;
          l2dot_temp_ptr1 += VECTOR_SIZE;

          temp_level_ptr += BLOCK_LENGTH;
          temp_level_int_ptr += BLOCK_LENGTH;
          temp_index_ptr += BLOCK_LENGTH;

          i_idx += BLOCK_LENGTH;
        }

        gradient_temp_ptr1 = gradient_temp_ptr;

        for (size_t d_outer = 0; d_outer < max_dims; d_outer++) {  // D + 2
          __m256d mm_element = _mm256_load_pd(alpha_padded_temp_ptr_ + j);
          __m256d mm_temp = _mm256_load_pd(gradient_temp_ptr1);
          mm_element = _mm256_mul_pd(mm_element, mm_temp);

          __m256d mm_element2 = _mm256_load_pd(alpha_padded_temp_ptr_ + j + VECTOR_SIZE);

          __m256d mm_temp2 = _mm256_load_pd(gradient_temp_ptr1 + temp_cols);

          mm_element2 = _mm256_mul_pd(mm_element2, mm_temp2);

          l2dot_temp_ptr1 = l2dot_temp_ptr;

          for (size_t d_inner = 0; d_inner < d_outer; d_inner++) {
            // element *= ((l2dot(i, j, d_inner)*(d_outer != d_inner)) + (gradient(i, j,
            // d_inner)*(d_outer == d_inner)));
            __m256d mm_temp = _mm256_load_pd(l2dot_temp_ptr1);
            mm_element = _mm256_mul_pd(mm_element, mm_temp);
            __m256d mm_temp2 = _mm256_load_pd(l2dot_temp_ptr1 + temp_cols);
            mm_element2 = _mm256_mul_pd(mm_element2, mm_temp2);

            l2dot_temp_ptr1 += VECTOR_SIZE;
          }

          l2dot_temp_ptr1 += VECTOR_SIZE;

          for (size_t d_inner = d_outer + 1; d_inner < max_dims; d_inner++) {
            // element *= ((l2dot(i, j, d_inner)*(d_outer != d_inner)) + (gradient(i, j,
            // d_inner)*(d_outer == d_inner)));
            __m256d mm_temp = _mm256_load_pd(l2dot_temp_ptr1);
            mm_element = _mm256_mul_pd(mm_element, mm_temp);
            __m256d mm_temp2 = _mm256_load_pd(l2dot_temp_ptr1 + temp_cols);
            mm_element2 = _mm256_mul_pd(mm_element2, mm_temp2);

            l2dot_temp_ptr1 += VECTOR_SIZE;
          }

          __m256d mm_lambda = _mm256_broadcast_sd(lambda_ptr_ + d_outer);
          mm_element = _mm256_mul_pd(mm_element, mm_lambda);
          mm_element2 = _mm256_mul_pd(mm_element2, mm_lambda);

          mm_result = _mm256_add_pd(mm_result, mm_element);
          mm_result2 = _mm256_add_pd(mm_result2, mm_element2);

          gradient_temp_ptr1 += VECTOR_SIZE;
        }

#if defined(STORE_MATRIX)
        double* operation_result_dest_ptr =
            operation_result_matrix_->getPointer() +
            (i - process_i_start) * operation_result_matrix_->getNcols();
        _mm256_store_pd(operation_result_dest_ptr + j, mm_result);
        _mm256_store_pd(operation_result_dest_ptr + j + VECTOR_SIZE, mm_result2);
#endif
      }

#if !defined(STORE_MATRIX)
      mm_result = _mm256_add_pd(mm_result, mm_result2);
      mm_result = _mm256_hadd_pd(mm_result, _mm256_setzero_pd());

      double s_result;
      __m256d hsum = _mm256_add_pd(mm_result, _mm256_permute2f128_pd(mm_result, mm_result, 0x1));
      _mm_store_sd(&s_result,
                   _mm_hadd_pd(_mm256_castpd256_pd128(hsum), _mm256_castpd256_pd128(hsum)));

      result_ptr_[i] += s_result;
#endif
    }
  }

#elif defined(__SSE4_2__) && !defined(__AVX__)

#pragma omp parallel
  {
    std::size_t padded_size = this->level_->getNcols();
    double* constants = this->constants_->getPointer();  // {0, 0.5, 2.0 / 3.0, 1, 2};
#if !defined(STORE_MATRIX)
    double* result_ptr_ = result.getPointer();
#endif
    double* level_ptr_ = this->level_->getPointer();
    double* level_int_ptr_ = this->level_int_->getPointer();
    double* index_ptr_ = this->index_->getPointer();

    double* gradient_temp_ptr = gradient_temp[omp_get_thread_num()]->getPointer();
    double* l2dot_temp_ptr = l2dot_temp[omp_get_thread_num()]->getPointer();
    double* lambda_ptr_ = this->lambda_->getPointer();
    size_t temp_cols = this->storage->getDimension() * VECTOR_SIZE;

    double* lcl_q_temp_ptr_ = lcl_q_->getPointer();
    double* lcl_q_inv_temp_ptr_ = lcl_q_inv_->getPointer();
    double* alpha_padded_temp_ptr_ = alpha_padded_->getPointer();

    size_t max_dims = this->storage->getDimension();
    size_t page_cap_rounded = max_dims * BLOCK_LENGTH;

    __m128d mm_half = _mm_loaddup_pd(constants + 1);
    __m128d mm_two_thirds = _mm_loaddup_pd(constants + 2);
    __m128d mm_one = _mm_loaddup_pd(constants + 3);
    __m128d mm_two = _mm_loaddup_pd(constants + 4);
    __m128d mm_abs = _mm_loaddup_pd(constants + 5);

    size_t thr_start;
    size_t thr_end;
    SGPP::parallel::PartitioningTool::getOpenMPPartitionSegment(process_i_start, process_i_end,
                                                                &thr_start, &thr_end);

    for (size_t i = thr_start; i < thr_end; i++) {
      __m128d mm_result = _mm_setzero_pd();
      __m128d mm_result2 = _mm_setzero_pd();

      size_t i_page = i / BLOCK_LENGTH;
      size_t i_offset = i_page * page_cap_rounded + i % BLOCK_LENGTH;

      double* temp_level_ptr = level_ptr_;
      double* temp_level_int_ptr = level_int_ptr_;
      double* temp_index_ptr = index_ptr_;

      for (size_t j = 0; j < padded_size; j += VECTOR_SIZE * REG_BCOUNT) {
#if defined(STORE_MATRIX)
        mm_result = _mm_setzero_pd();
        mm_result2 = _mm_setzero_pd();
#endif

        double* gradient_temp_ptr1 = gradient_temp_ptr;
        double* l2dot_temp_ptr1 = l2dot_temp_ptr;

        size_t i_idx = i_offset;

        for (size_t dim = 0; dim < max_dims; dim++) {
          __m128d mm_lcl_q = _mm_loaddup_pd(lcl_q_temp_ptr_ + dim);
          __m128d mm_lcl_q_inv = _mm_loaddup_pd(lcl_q_inv_temp_ptr_ + dim);

          __m128d mm_lid = _mm_loaddup_pd(level_ptr_ + i_idx);
          __m128d mm_iid = _mm_loaddup_pd(index_ptr_ + i_idx);
          __m128d mm_ljd = _mm_load_pd(temp_level_ptr);
          __m128d mm_ijd = _mm_load_pd(temp_index_ptr);

          __m128d mm_doGrad =
              _mm_and_pd(_mm_cmpeq_pd(mm_lid, mm_ljd), _mm_cmpeq_pd(mm_iid, mm_ijd));  // +3

          __m128d mm_grad = _mm_mul_pd(mm_lid, _mm_and_pd(mm_two, mm_doGrad));  // 1+1

          mm_grad = _mm_mul_pd(mm_grad, mm_lcl_q_inv);  // 1

          _mm_store_pd(gradient_temp_ptr1, mm_grad);

          __m128d mm_in_lid = _mm_loaddup_pd(level_int_ptr_ + i_idx);
          __m128d mm_in_ljd = _mm_load_pd(temp_level_int_ptr);

          __m128d mm_res_one = _mm_mul_pd(
              mm_two_thirds, _mm_and_pd(mm_in_lid, _mm_cmpeq_pd(mm_iid, mm_ijd)));  // 1+2

          __m128d mm_selector = _mm_cmple_pd(mm_lid, mm_ljd);  // +6
          __m128d mm_i1d = _mm_blendv_pd(mm_iid, mm_ijd, mm_selector);
          __m128d mm_in_l1d = _mm_blendv_pd(mm_in_lid, mm_in_ljd, mm_selector);
          __m128d mm_in_l2d = _mm_blendv_pd(mm_in_ljd, mm_in_lid, mm_selector);
          __m128d mm_i2d = _mm_blendv_pd(mm_ijd, mm_iid, mm_selector);
          __m128d mm_l2d = _mm_blendv_pd(mm_ljd, mm_lid, mm_selector);

          __m128d mm_q = _mm_mul_pd(_mm_sub_pd(mm_i1d, mm_one), mm_in_l1d);  // 2 flop
          __m128d mm_p = _mm_mul_pd(_mm_add_pd(mm_i1d, mm_one), mm_in_l1d);  // 2 flop
          __m128d mm_overlap = _mm_cmplt_pd(
              _mm_max_pd(mm_q, _mm_mul_pd(_mm_sub_pd(mm_i2d, mm_one), mm_in_l2d)),
              _mm_min_pd(mm_p, _mm_mul_pd(_mm_add_pd(mm_i2d, mm_one), mm_in_l2d)));  // 6+1

          __m128d mm_temp_res = _mm_sub_pd(
              _mm_sub_pd(mm_two,
                         _mm_and_pd(mm_abs, (_mm_sub_pd(_mm_mul_pd(mm_l2d, mm_q), mm_i2d)))),
              _mm_and_pd(mm_abs, (_mm_sub_pd(_mm_mul_pd(mm_l2d, mm_p), mm_i2d))));  // 8 flops
          mm_temp_res = _mm_mul_pd(mm_temp_res, _mm_mul_pd(mm_half, mm_in_l1d));    // 2 flops
          __m128d mm_res_two = _mm_and_pd(mm_temp_res, mm_overlap);  // Now mask result //+1
          mm_selector = _mm_cmpneq_pd(mm_lid, mm_ljd);               // +1

          __m128d mm_val = _mm_blendv_pd(mm_res_one, mm_res_two, mm_selector);  // +1
          mm_val = _mm_mul_pd(mm_val, mm_lcl_q);                                // 1 flop
          _mm_store_pd(l2dot_temp_ptr1, mm_val);

          ////////////////////////////////////////////////////////
          __m128d mm_ljd2 = _mm_load_pd(temp_level_ptr + VECTOR_SIZE);
          __m128d mm_ijd2 = _mm_load_pd(temp_index_ptr + VECTOR_SIZE);
          __m128d mm_doGrad2 =
              _mm_and_pd(_mm_cmpeq_pd(mm_lid, mm_ljd2), _mm_cmpeq_pd(mm_iid, mm_ijd2));  // 1 // +2

          __m128d mm_grad2 = _mm_mul_pd(mm_lid, _mm_and_pd(mm_two, mm_doGrad2));  // 2

          mm_grad2 = _mm_mul_pd(mm_grad2, mm_lcl_q_inv);  // 1
          _mm_store_pd(gradient_temp_ptr1 + temp_cols, mm_grad2);

          __m128d mm_in_ljd2 = _mm_load_pd(temp_level_int_ptr + VECTOR_SIZE);

          __m128d mm_res_one2 = _mm_mul_pd(
              mm_two_thirds, _mm_and_pd(mm_in_lid, _mm_cmpeq_pd(mm_iid, mm_ijd2)));  // 2 // +1
          __m128d mm_selector2 = _mm_cmple_pd(mm_lid, mm_ljd2);
          __m128d mm_i1d2 = _mm_blendv_pd(mm_iid, mm_ijd2, mm_selector2);
          __m128d mm_in_l1d2 = _mm_blendv_pd(mm_in_lid, mm_in_ljd2, mm_selector2);
          __m128d mm_i2d2 = _mm_blendv_pd(mm_ijd2, mm_iid, mm_selector2);
          __m128d mm_l2d2 = _mm_blendv_pd(mm_ljd2, mm_lid, mm_selector2);
          __m128d mm_in_l2d2 = _mm_blendv_pd(mm_in_ljd2, mm_in_lid, mm_selector2);

          __m128d mm_q2 = _mm_mul_pd(_mm_sub_pd(mm_i1d2, mm_one), mm_in_l1d2);  // 2 flop
          __m128d mm_p2 = _mm_mul_pd(_mm_add_pd(mm_i1d2, mm_one), mm_in_l1d2);  // 2 flop

          __m128d mm_overlap2 =
              _mm_cmplt_pd(_mm_max_pd(mm_q2, _mm_mul_pd(_mm_sub_pd(mm_i2d2, mm_one), mm_in_l2d2)),
                           _mm_min_pd(mm_p2, _mm_mul_pd(_mm_add_pd(mm_i2d2, mm_one),
                                                        mm_in_l2d2)));  // 6 flop // +1

          __m128d mm_temp_res2 = _mm_sub_pd(
              _mm_sub_pd(mm_two,
                         _mm_and_pd(mm_abs, (_mm_sub_pd(_mm_mul_pd(mm_l2d2, mm_q2), mm_i2d2)))),
              _mm_and_pd(mm_abs, (_mm_sub_pd(_mm_mul_pd(mm_l2d2, mm_p2),
                                             mm_i2d2))));  // 8 flops

          mm_temp_res2 = _mm_mul_pd(mm_temp_res2, _mm_mul_pd(mm_half, mm_in_l1d2));  // 2 flops
          __m128d mm_res_two2 = _mm_and_pd(mm_temp_res2, mm_overlap2);  // Now mask result //1 flop

          mm_selector2 = _mm_cmpneq_pd(mm_lid, mm_ljd2);

          __m128d mm_val2 = _mm_blendv_pd(mm_res_one2, mm_res_two2, mm_selector2);

          mm_val2 = _mm_mul_pd(mm_val2, mm_lcl_q);  // 1 flop

          _mm_store_pd(l2dot_temp_ptr1 + temp_cols, mm_val2);

          gradient_temp_ptr1 += VECTOR_SIZE;
          l2dot_temp_ptr1 += VECTOR_SIZE;

          temp_level_ptr += BLOCK_LENGTH;
          temp_level_int_ptr += BLOCK_LENGTH;
          temp_index_ptr += BLOCK_LENGTH;

          i_idx += BLOCK_LENGTH;
        }

        gradient_temp_ptr1 = gradient_temp_ptr;

        for (size_t d_outer = 0; d_outer < max_dims; d_outer++) {  // D + 2
          __m128d mm_element = _mm_load_pd(alpha_padded_temp_ptr_ + j);
          __m128d mm_temp = _mm_load_pd(gradient_temp_ptr1);
          mm_element = _mm_mul_pd(mm_element, mm_temp);

          __m128d mm_element2 = _mm_load_pd(alpha_padded_temp_ptr_ + j + VECTOR_SIZE);

          __m128d mm_temp2 = _mm_load_pd(gradient_temp_ptr1 + temp_cols);

          mm_element2 = _mm_mul_pd(mm_element2, mm_temp2);

          l2dot_temp_ptr1 = l2dot_temp_ptr;

          for (size_t d_inner = 0; d_inner < d_outer; d_inner++) {
            // element *= ((l2dot(i, j, d_inner)*(d_outer != d_inner)) + (gradient(i, j,
            // d_inner)*(d_outer == d_inner)));
            __m128d mm_temp = _mm_load_pd(l2dot_temp_ptr1);
            mm_element = _mm_mul_pd(mm_element, mm_temp);
            __m128d mm_temp2 = _mm_load_pd(l2dot_temp_ptr1 + temp_cols);
            mm_element2 = _mm_mul_pd(mm_element2, mm_temp2);

            l2dot_temp_ptr1 += VECTOR_SIZE;
          }

          l2dot_temp_ptr1 += VECTOR_SIZE;

          for (size_t d_inner = d_outer + 1; d_inner < max_dims; d_inner++) {
            // element *= ((l2dot(i, j, d_inner)*(d_outer != d_inner)) + (gradient(i, j,
            // d_inner)*(d_outer == d_inner)));
            __m128d mm_temp = _mm_load_pd(l2dot_temp_ptr1);
            mm_element = _mm_mul_pd(mm_element, mm_temp);
            __m128d mm_temp2 = _mm_load_pd(l2dot_temp_ptr1 + temp_cols);
            mm_element2 = _mm_mul_pd(mm_element2, mm_temp2);

            l2dot_temp_ptr1 += VECTOR_SIZE;
          }

          __m128d mm_lambda = _mm_loaddup_pd(lambda_ptr_ + d_outer);
          mm_element = _mm_mul_pd(mm_element, mm_lambda);
          mm_element2 = _mm_mul_pd(mm_element2, mm_lambda);

          mm_result = _mm_add_pd(mm_result, mm_element);
          mm_result2 = _mm_add_pd(mm_result2, mm_element2);

          gradient_temp_ptr1 += VECTOR_SIZE;
        }

#if defined(STORE_MATRIX)
        double* operation_result_dest_ptr =
            operation_result_matrix_->getPointer() +
            (i - process_i_start) * operation_result_matrix_->getNcols();
        _mm_store_pd(operation_result_dest_ptr + j, mm_result);
        _mm_store_pd(operation_result_dest_ptr + j + VECTOR_SIZE, mm_result2);
#endif
      }

#if !defined(STORE_MATRIX)
      mm_result = _mm_add_pd(mm_result, mm_result2);
      mm_result = _mm_hadd_pd(mm_result, _mm_setzero_pd());

      double s_result = 0.0;
      _mm_store_sd(&s_result, mm_result);

      result_ptr_[i] += s_result;
#endif
    }
  }

#else
#pragma omp parallel
  {
    double* gradient_temp_ptr = gradient_temp[omp_get_thread_num()]->getPointer();
    double* l2dot_temp_ptr = l2dot_temp[omp_get_thread_num()]->getPointer();

    size_t thr_start;
    size_t thr_end;
    SGPP::parallel::PartitioningTool::getOpenMPPartitionSegment(process_i_start, process_i_end,
                                                                &thr_start, &thr_end);

    for (size_t i = thr_start; i < thr_end; i++) {
      for (size_t j = 0; j < this->storage->getSize(); j++) {
        double temp = 0.0;

        for (size_t d = 0; d < this->storage->getDimension(); d++) {
          l2dot_temp_ptr[d] = l2dot(i, j, d);
          gradient_temp_ptr[d] = gradient(i, j, d);
        }

        for (size_t d_outer = 0; d_outer < this->storage->getDimension(); d_outer++) {
          double element = alpha[j];

          for (size_t d_inner = 0; d_inner < this->storage->getDimension(); d_inner++) {
            element *= ((l2dot_temp_ptr[d_inner] * (d_outer != d_inner)) +
                        (gradient_temp_ptr[d_inner] * (d_outer == d_inner)));
          }

          temp += (this->lambda_->get(d_outer) * element);
        }

        result[i] += temp;

#if defined(STORE_MATRIX)
        double* operation_result_dest_ptr =
            operation_result_matrix_->getPointer() +
            (i - process_i_start) * operation_result_matrix_->getNcols();
        operation_result_dest_ptr[j] = temp;
#endif
      }
    }
  }
#endif

#if defined(STORE_MATRIX)
  }

  std::size_t original_size = alpha.getSize();
  memcpy(alpha_padded_->getPointer(), alpha.getPointer(), original_size * sizeof(double));

#pragma omp parallel
  {
    size_t thr_start;
    size_t thr_end;
    SGPP::parallel::PartitioningTool::getOpenMPPartitionSegment(process_i_start, process_i_end,
                                                                &thr_start, &thr_end);

    for (size_t i = thr_start; i < thr_end; i++) {
      double* operation_result_dest_ptr =
          operation_result_matrix_->getPointer() +
          (i - process_i_start) * operation_result_matrix_->getNcols();

      double element = 0.0;

#if defined(__MIC__)
#pragma prefetch
#endif

      for (size_t j = 0; j < storage->getSize(); ++j) {
        element += alpha[j] * *(operation_result_dest_ptr + j);
      }

      result[i] = element;
    }
  }
#endif

#ifdef USE_MPI
  double* result_ptr = result.getPointer();
  /*
  MPI_Alltoallv(result_ptr, send_size.data(), send_start.data(), MPI_DOUBLE,
           result_ptr, recv_size.data(), recv_start.data(), MPI_DOUBLE,
           MPI_COMM_WORLD);
  */
  MPI_Allgatherv(MPI_IN_PLACE, send_size[0], MPI_DOUBLE, result_ptr, recv_size.data(),
                 recv_start.data(), MPI_DOUBLE, MPI_COMM_WORLD);
/*
MPI_Allreduce(MPI_IN_PLACE, result_ptr, (int)result.getSize(), MPI_DOUBLE,
           MPI_SUM, MPI_COMM_WORLD);
*/
#endif

  all_time += stopWatch.stop();
  all_iterations += 1.0;
}
}  // namespace parallel
}  // namespace SGPP
