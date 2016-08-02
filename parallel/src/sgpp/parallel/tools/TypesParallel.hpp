// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TYPESPARALLEL_HPP
#define TYPESPARALLEL_HPP

#include <sgpp/globaldef.hpp>

namespace sgpp {

namespace parallel {

enum VectorizationType {
  X86SIMD,
  OpenCL,
  Hybrid_X86SIMD_OpenCL,
  MIC,
  Hybrid_X86SIMD_MIC,
  CUDA
};

enum MPIType {
  MPINone,
  MPIAllreduce,
  MPIAsync,
  MPIOnesided,
  MPIAlltoallv,
  MPITrueAsync,
  MPIBigdata
};
}  // namespace parallel
}  // namespace sgpp

#endif /* TYPESPARALLEL_HPP */
