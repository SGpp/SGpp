// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef KERNELBASE_HPP
#define KERNELBASE_HPP

#include <sgpp/base/grid/GridStorage.hpp>

// #define CHECK_KERNEL_CALLS
#ifdef CHECK_KERNEL_CALLS
#include <sgpp/base/exception/operation_exception.hpp>

// helper macros
#define ASSERT_EQUAL(arg1, arg2)                                                             \
  {                                                                                          \
    if ((arg1) != (arg2)) {                                                                  \
      std::cerr << #arg1 << " and " << #arg2 << " are not equal: " << arg1 << " != " << arg2 \
                << " (file:" << __FILE__ << ", line:" << __LINE__ << ")" << std::endl;       \
      throw sgpp::base::operation_exception("values " #arg1 " and " #arg2 " are not equal"); \
    }                                                                                        \
  }
#define ASSERT_ALIGNMENT(arg, alignment)                                                       \
  {                                                                                            \
    if (((arg) % (alignment)) != 0) {                                                          \
      std::cout << #arg << "(" << arg << ") not aligned to " << #alignment << "(" << alignment \
                << ") (file:" << __FILE__ << ", line:" << __LINE__ << ")" << std::endl;        \
      throw sgpp::base::operation_exception("argument " #arg " not aligned!");                 \
    }                                                                                          \
  }
#define ASSERT_GEQ_THAN(arg1, arg2)                                                     \
  {                                                                                     \
    if ((arg1) < (arg2)) {                                                              \
      std::cout << #arg1 << "(" << arg1 << ") is smaller than " << #arg2 << "(" << arg2 \
                << ") (file:" << __FILE__ << ", line:" << __LINE__ << ")" << std::endl; \
      throw sgpp::base::operation_exception(#arg1 " is smaller than " #arg2);           \
    }                                                                                   \
  }
#define ASSERT_LEQ_THAN(arg1, arg2)                                                     \
  {                                                                                     \
    if ((arg1) > (arg2)) {                                                              \
      std::cout << #arg1 << "(" << arg1 << ") is greater than " << #arg2 << "(" << arg2 \
                << ") (file:" << __FILE__ << ", line:" << __LINE__ << ")" << std::endl; \
      throw sgpp::base::operation_exception(#arg1 " is greater than " #arg2);           \
    }                                                                                   \
  }

#define CHECK_INDEX_ARG(arg, min, max) \
  {                                    \
    ASSERT_GEQ_THAN(arg, min)          \
    ASSERT_LEQ_THAN(arg, max)          \
  }
#define CHECK_DATASET_AND_SOURCE(dataset, source)                   \
  {                                                                 \
    /*ASSERT_ALIGNMENT(dataset->getNcols(), getChunkDataPoints())*/ \
    ASSERT_EQUAL(dataset->getNcols(), source.getSize())             \
  }
#define CHECK_DATASET_AND_RESULT(dataset, result)                   \
  {                                                                 \
    /*ASSERT_ALIGNMENT(dataset->getNcols(), getChunkDataPoints())*/ \
    ASSERT_EQUAL(dataset->getNcols(), result.getSize())             \
  }

// use only the following two
#define CHECK_ARGS_MULT(level, dataset, result, s_grid, e_grid, s_data, e_data) \
  {                                                                             \
    CHECK_INDEX_ARG(s_grid, 0, level->getNrows());                              \
    CHECK_INDEX_ARG(e_grid, 0, level->getNrows());                              \
    CHECK_INDEX_ARG(s_data, 0, result.getSize());                               \
    CHECK_INDEX_ARG(e_data, 0, result.getSize());                               \
    ASSERT_ALIGNMENT(e_data - s_data, getChunkDataPoints());                    \
    CHECK_DATASET_AND_RESULT(dataset, result);                                  \
  }
#define CHECK_ARGS_MULTTRANSPOSE(level, dataset, source, s_grid, e_grid, s_data, e_data) \
  {                                                                                      \
    CHECK_INDEX_ARG(s_grid, 0, level->getNrows());                                       \
    CHECK_INDEX_ARG(e_grid, 0, level->getNrows());                                       \
    CHECK_INDEX_ARG(s_data, 0, source.getSize());                                        \
    CHECK_INDEX_ARG(e_data, 0, source.getSize());                                        \
    ASSERT_ALIGNMENT(e_data - s_data, getChunkDataPoints());                             \
    CHECK_DATASET_AND_SOURCE(dataset, source);                                           \
  }
#else
#define CHECK_ARGS_MULT(level, dataset, result, s_grid, e_grid, s_data, e_data)
#define CHECK_ARGS_MULTTRANSPOSE(level, dataset, source, s_grid, e_grid, s_data, e_data)
#endif

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {
enum KernelType { Standard, Mask };

}  // namespace parallel
}  // namespace sgpp

#endif  // KERNELBASE_HPP
