// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/OperationMultiEvalCuda/MortonOrder.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <stdint.h>

#include <sgpp/datadriven/operation/hash/OperationMultiEvalCuda/MortonOrderKernel.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultiEvalCuda/basicCuda.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultiEvalCuda/cudaHelper.hpp>

///@cond DOXY_IGNORE // NOLINT()
namespace sgpp {
namespace datadriven {
namespace OpMultiEvalCudaDetail {

/// Constructor. Generates the identity
MortonOrder::MortonOrder(size_t size) {
  permutation.resize(size);
  for (size_t i = 0; i < size; ++i) {
    permutation[i] = i;
  }
}
/// Constructor. Generates the permuation list on the GPU
MortonOrder::MortonOrder(sgpp::base::DataMatrix& matrix) {
  size_t DIM = matrix.getNcols();
  HostDevPtr<size_t> perm;
  HostDevPtr<double> pos;
  size_t N = matrix.getNrows();
  perm.HostAlloc(N);
  pos.HostAlloc(N * DIM);
  // Init permuation list without any permuation
  for (size_t i = 0; i < N; ++i) {
    for (size_t d = 0; d < DIM; ++d) {
      perm[i] = i;
      pos[i * DIM + d] = matrix.get(i, d);
    }
  }
  perm.DevAlloc();
  pos.DevAlloc();
  perm.CopyToDev();
  pos.CopyToDev();
  CudaCheckError();

  // Compute permuation
  zorder(pos.dev(), perm.dev(), N, DIM);
  CudaCheckError();

  // Copy back to host
  perm.CopyToHost();
  permutation.resize(N);
  for (size_t i = 0; i < N; ++i) {
    permutation[i] = perm[i];
  }
}

/// Apply permuation list
void MortonOrder::orderDataMatrix(sgpp::base::DataMatrix& matrix) const {
  if (matrix.getNrows() != permutation.size())
    throw sgpp::base::operation_exception("dimensions do not match!");
  sgpp::base::DataMatrix origin(matrix);
  for (size_t i = 0; i < matrix.getNrows(); ++i) {
    for (size_t d = 0; d < matrix.getNcols(); ++d) {
      matrix(i, d) = origin(permutation[i], d);
    }
  }
}
/// Apply permuation list
void MortonOrder::orderDataMatrix(sgpp::base::DataMatrix& matrix, double* dest) const {
  if (matrix.getNrows() != permutation.size())
    throw sgpp::base::operation_exception("dimensions do not match!");
  for (size_t i = 0; i < matrix.getNrows(); ++i) {
    for (size_t d = 0; d < matrix.getNcols(); ++d) {
      dest[i * matrix.getNcols() + d] = matrix(permutation[i], d);
    }
  }
}
/// Apply permuation list
void MortonOrder::orderDataVector(sgpp::base::DataVector& data) const {
  if (data.getSize() != permutation.size())
    throw sgpp::base::operation_exception("dimensions do not match!");
  sgpp::base::DataVector origin(data);
  for (size_t i = 0; i < data.getSize(); ++i) {
    data[i] = origin[permutation[i]];
  }
}
/// Apply permuation list
void MortonOrder::orderDataVector(sgpp::base::DataVector& data, double* dest) const {
  if (data.getSize() != permutation.size())
    throw sgpp::base::operation_exception("dimensions do not match!");
  for (size_t i = 0; i < data.getSize(); ++i) {
    dest[i] = data[permutation[i]];
  }
}

/// Apply permuation list in reverse
void MortonOrder::restoreDataMatrix(sgpp::base::DataMatrix& matrix) const {
  if (matrix.getNrows() != permutation.size())
    throw sgpp::base::operation_exception("dimensions do not match!");
  sgpp::base::DataMatrix origin(matrix);
  for (size_t i = 0; i < matrix.getNrows(); ++i) {
    for (size_t d = 0; d < matrix.getNcols(); ++d) {
      matrix(permutation[i], d) = origin(i, d);
    }
  }
}
/// Apply permuation list in reverse
void MortonOrder::restoreDataMatrix(sgpp::base::DataMatrix& matrix, double* src) const {
  if (matrix.getNrows() != permutation.size())
    throw sgpp::base::operation_exception("dimensions do not match!");
  for (size_t i = 0; i < matrix.getNrows(); ++i) {
    for (size_t d = 0; d < matrix.getNcols(); ++d) {
      matrix(permutation[i], d) = src[i * matrix.getNcols() + d];
    }
  }
}
/// Apply permuation list in reverse
void MortonOrder::restoreDataVector(sgpp::base::DataVector& data) const {
  if (data.getSize() != permutation.size())
    throw sgpp::base::operation_exception("dimensions do not match!");
  sgpp::base::DataVector origin(data);
  for (size_t i = 0; i < data.getSize(); ++i) {
    data[permutation[i]] = origin[i];
  }
}
/// Apply permuation list in reverse
void MortonOrder::restoreDataVector(sgpp::base::DataVector& data, double* src) const {
  if (data.getSize() != permutation.size())
    throw sgpp::base::operation_exception("dimensions do not match!");
  for (size_t i = 0; i < data.getSize(); ++i) {
    data[permutation[i]] = src[i];
  }
}
/// Check if permutation is identity
bool MortonOrder::isIdentity() const {
  for (size_t i = 0; i < permutation.size(); ++i)
    if (permutation[i] != i) return false;
  return true;
}

}  // namespace OpMultiEvalCudaDetail
}  // namespace datadriven
}  // namespace sgpp
///@endcond // NOLINT()
