#include "OperationMultiEvalCuda.hpp"
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>

#include "cudaHelper.hpp"
#include "basicCuda.hpp"
#include "MultiEvalKernel.hpp"
#include "MortonOrderKernel.hpp"
#include "consts.hpp"
#include <stdint.h>

namespace sgpp {
namespace datadriven {

/// Construct the class
OperationMultiEvalCuda::OperationMultiEvalCuda(base::Grid& grid, base::DataMatrix& dataset,
    const uint32_t& grad) : base::OperationMultipleEval(grid, dataset), _ordered(false) {
  polygrad = grad;
  prepare();
}

/// Construct the class
OperationMultiEvalCuda::OperationMultiEvalCuda(base::Grid& grid, base::DataMatrix& dataset, 
    const uint32_t& grad, bool ordered) : base::OperationMultipleEval(grid, dataset), _ordered(ordered) {
  polygrad = grad;
  prepare();
}

/// Destroy everything ;)
OperationMultiEvalCuda::~OperationMultiEvalCuda() {
}

/// Standard evaluation
void OperationMultiEvalCuda::mult(sgpp::base::DataVector& source, sgpp::base::DataVector& result) {
  myTimer.start();
  // Copy stuff to device
  for (uint32_t i = 0; i < _N; ++i) {
    alpha[i] = source[i];
  }
  for (uint32_t i = _N; i < N; ++i) {
    alpha[i] = 0;
  }
  alpha.CopyToDev();
  
  // Call eval kernel
  evalCuda(data.dev(),alpha.dev(),node.dev(),pos.dev(),M,maxlevel,levellimit.dev(),_b[maxlevel][DIM],
    subs.dev());
    
  // Copy result back to host
  CudaCheckError();
  data.CopyToHost();
  for (uint32_t j = 0; j < M; ++j) {
    result[j] = data[j];
  }
  duration = myTimer.stop();
}

/// Transposed evaluation
void OperationMultiEvalCuda::multTranspose(sgpp::base::DataVector& source, 
    sgpp::base::DataVector& result) {
  myTimer.start();
  // Copy stuff to device
  for (uint32_t j = 0; j < M; ++j) {
    data[j] = source[j];
  }
  data.CopyToDev();
  
  // Call transposed eval kernel
  transposedCuda(alpha.dev(), node.dev(), pos.dev(), data.dev(), streamlimit.dev(), M, _N);
  
  // Copy result back to host
  CudaCheckError();
  alpha.CopyToHost();
  for (uint32_t i = 0; i < _N; ++i) {
    result[i] = alpha[i];
  }
  duration = myTimer.stop();
}

/// Transposed evaluation with additional FMA
void OperationMultiEvalCuda::multTransposeFMA(sgpp::base::DataVector& source,
    sgpp::base::DataVector& prev, double lambda, sgpp::base::DataVector& result) {
  myTimer.start();
  // Copy stuff to device
  HostDevPtr<double> b;
  b.HostAlloc(N);
  for (uint32_t j = 0; j < M; ++j) {
    data[j] = source[j];
  }
  for (uint32_t i = 0; i < _N; ++i) {
     b[i] = prev[i];
  }
  for (uint32_t i = _N; i < N; ++i) {
     b[i] = 0;
  }
  b.DevAlloc();
  data.CopyToDev();
  b.CopyToDev();
  CudaCheckError();
  
  // Call transposed evaluation kernel
  transposedCuda(alpha.dev(), node.dev(), pos.dev(), data.dev(), streamlimit.dev(), b.dev(),
    lambda, M, _N);
    
  // Copy result back to host
  CudaCheckError();
  alpha.CopyToHost();
  for (uint32_t i = 0; i < _N; ++i) {
    result[i] = alpha[i];
  }
  duration = myTimer.stop();
}

/// Preprocessing
void OperationMultiEvalCuda::prepare() {
  // Prepare grid
  sgpp::base::GridStorage gs(grid.getStorage());
  DIM = gs.getDimension();
  _N = gs.getSize();
  M = dataset.getNrows();
  N = (_N + CUDA_BLOCKSIZE - 1)&(0xFFFFFFFF - CUDA_BLOCKSIZE + 1);
  
  node.HostAlloc(DIM*N);
  alpha.HostAlloc(N);
  streamlimit.HostAlloc(N);
  
  // Copy SG++ grid structure to graph based datastructure on the GPU
  sgpp::base::GridPoint gp;
  maxlevel = 1;
  for (uint32_t i = 0; i < _N; ++i) {
    if (!_ordered) {
      streamlimit[i].lower = 0;
      streamlimit[i].upper = M;
    }
    uint32_t ml = 1;
    for (uint32_t d = 0; d < DIM; ++d) {
      gp = gs.getPoint(i);
      ml += gp.getLevel(d);
      node[i*DIM + d].level = gp.getLevel(d);
      node[i*DIM + d].index = gp.getIndex(d);
      node[i*DIM + d].grad = (gp.getLevel(d)>=polygrad?polygrad - 1:gp.getLevel(d));
      node[i*DIM + d].level2 = double(1<<gp.getLevel(d));
      node[i*DIM + d].x = double(gp.getIndex(d))/node[i*DIM + d].level2;
      gp.getLeftChild(d);
      if (gs.isContaining(gp))
        node[i*DIM + d].child[0] = gs.getSequenceNumber(gp);
      else
        node[i*DIM + d].child[0] = GRID_END;
      gp = gs.getPoint(i);
      gp.getRightChild(d);
      if (gs.isContaining(gp))
        node[i*DIM + d].child[1] = gs.getSequenceNumber(gp);
      else
        node[i*DIM + d].child[1] = GRID_END;
    }
    ml -= DIM;
    if (ml > maxlevel) maxlevel = ml;
  }
  
  // Fill up to multiple of block size
  for (uint32_t i = _N; i < N; ++i) {
    if (!_ordered) {
      streamlimit[i].lower = 0;
      streamlimit[i].upper = M;
    }
    for (uint32_t d = 0; d < DIM; ++d) {
      node[i*DIM + d].level = node[(_N - 1)*DIM + d].level;
      node[i*DIM + d].index = node[(_N - 1)*DIM + d].index;
      node[i*DIM + d].grad = node[(_N - 1)*DIM + d].grad;
      node[i*DIM + d].level2 = node[(_N - 1)*DIM + d].level2;
      node[i*DIM + d].x = node[(_N - 1)*DIM + d].x;
      node[i*DIM + d].child[0] = node[(_N - 1)*DIM + d].child[0];
      node[i*DIM + d].child[1] = node[(_N - 1)*DIM + d].child[1];
    }
  }
  
  // Prepare subspace index
  // Computes all subspaces with levelsum = maxlevel + d - 1 
  subs.HostAlloc(_b[maxlevel][DIM]*DIM);
  uint32_t gnode[DIM];
	uint32_t sum = 1;
	uint32_t cnt = 0;
	for (uint32_t d=0;d<DIM;++d) gnode[d] = 1;
	for (uint32_t s=0;s<_b[maxlevel][DIM+1];++s) {
		if (sum == maxlevel) {
      for (uint32_t d=0;d<DIM;++d)
        subs[cnt*DIM + d] = gnode[d];
      cnt++;
    }
    for (uint32_t d=0;d<DIM;++d) {
      ++gnode[d];
      ++sum;
      if (sum <= maxlevel) break;
      sum -= gnode[d] - 1;
      gnode[d] = 1;
    }
	}
  
  // Prepare data
  pos.HostAlloc(DIM*M);
  data.HostAlloc(M);
  levellimit.HostAlloc(DIM*M); 
 
  for (uint32_t j = 0; j < M; ++j) {
    for (uint32_t d = 0; d < DIM; ++d) {
      pos[j*DIM + d] = dataset.get(j,d);
    }
  }
  
  // Preprocessing
  
  // Alloc all the memory on the device
  node.DevAlloc();
  alpha.DevAlloc();
  pos.DevAlloc();
  data.DevAlloc();
  streamlimit.DevAlloc();
  levellimit.DevAlloc();
  subs.DevAlloc();
  
  node.CopyToDev();
  pos.CopyToDev();
  subs.CopyToDev();
  
  // If dataset is not ordered along Morton curve
  if (!_ordered) {
    streamlimit.CopyToDev();
  }
  
  // Compute subspace limits
  preprocessCuda(node.dev(), pos.dev(), levellimit.dev(), maxlevel, M, DIM);
  CudaCheckError();
  
  // If dataset is aligned along Morton curve
  if (_ordered) {
    // Compute stream boundary limitations
    streamboundCuda(pos.dev(), node.dev(), streamlimit.dev(), M, N);
    CudaCheckError();
  }
  
  streamlimit.CopyToHost();
  levellimit.CopyToHost();
  
  isPrepared = true;
}

/// Returns last measured time
double OperationMultiEvalCuda::getDuration() {
  return duration;
}

}  // datadriven
}  // sgpp
