#include "MortonOrder.hpp"
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <stdint.h>

#include "cudaHelper.hpp"
#include "basicCuda.hpp"
#include "MortonOrderKernel.hpp"

namespace sgpp {
namespace datadriven {

/// Constructor. Generates the permuation list on the GPU
MortonOrder::MortonOrder (sgpp::base::DataMatrix& matrix) {
  uint32_t DIM = matrix.getNcols();
  HostDevPtr<size_t> perm;
  HostDevPtr<double> pos;
  perm.HostAlloc(matrix.getNrows());
  pos.HostAlloc(matrix.getSize());
  // Init permuation list without any permuation
  for (uint32_t i=0;i<matrix.getNrows();++i) {
    for (uint32_t d=0;d<DIM;++d) {
      perm[i] = i;
      pos[i*DIM + d] = matrix.get(i,d);
    }
  }
  perm.DevAlloc();
  pos.DevAlloc();
  perm.CopyToDev();
  pos.CopyToDev();
  CudaCheckError();
  
  // Compute permuation
  zorder(pos.dev(), perm.dev(), matrix.getNrows(), DIM);
  CudaCheckError();
  
  // Copy back to host
  perm.CopyToHost();
  permutation.resize(matrix.getNrows());
  for (uint32_t i=0;i<matrix.getNrows();++i) {
    permutation[i] = perm[i];
  }
}

/// Apply permuation list
void MortonOrder::orderDataMatrix(sgpp::base::DataMatrix& matrix) const {
  if (matrix.getNrows() != permutation.size())
    throw sgpp::base::operation_exception("dimensions do not match!");
  sgpp::base::DataMatrix origin(matrix);
  for (uint32_t i=0;i<matrix.getNrows();++i) {
    for (uint32_t d=0;d<matrix.getNcols();++d) {
      matrix(i,d) = origin(permutation[i],d);
    }
  }
}
/// Apply permuation list
void MortonOrder::orderDataVector(sgpp::base::DataVector& data) const {
  if (data.getSize() != permutation.size())
    throw sgpp::base::operation_exception("dimensions do not match!");
  sgpp::base::DataVector origin(data);
  for (uint32_t i=0;i<data.getSize();++i) {
    data[i] = origin[permutation[i]];
  }
}

/// Apply permuation list in reverse
void MortonOrder::restoreDataMatrix(sgpp::base::DataMatrix& matrix) const {
  if (matrix.getNrows() != permutation.size())
    throw sgpp::base::operation_exception("dimensions do not match!");
  sgpp::base::DataMatrix origin(matrix);
  for (uint32_t i=0;i<matrix.getNrows();++i) {
    for (uint32_t d=0;d<matrix.getNcols();++d) {
      matrix(permutation[i],d) = origin(i,d);
    }
  }
}
/// Apply permuation list in reverse
void MortonOrder::restoreDataVector(sgpp::base::DataVector& data) const {
  if (data.getSize() != permutation.size())
    throw sgpp::base::operation_exception("dimensions do not match!");
  sgpp::base::DataVector origin(data);
  for (uint32_t i=0;i<data.getSize();++i) {
    data[permutation[i]] = origin[i];
  }
}

}  // datadriven
}  // sgpp
