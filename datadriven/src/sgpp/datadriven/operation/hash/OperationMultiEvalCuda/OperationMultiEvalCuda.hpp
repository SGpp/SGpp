// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIEVALCUDA_HPP
#define OPERATIONMULTIEVALCUDA_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <stdint.h>

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/operation/hash/OperationMultiEvalCuda/MortonOrder.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultiEvalCuda/basicCuda.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultiEvalCuda/consts.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultiEvalCuda/cudaHelper.hpp>

namespace sgpp {
namespace datadriven {

/** OperationMultipleEval for polynomial basis functions (grad >= 2) using CUDA on grids
 *  without boundary nodes
 */
class OperationMultiEvalCuda : public base::OperationMultipleEval {
 private:
  bool _ordered;

 protected:
  OpMultiEvalCudaDetail::MortonOrder* zorder;

  base::SGppStopwatch myTimer;
  double duration;

  OpMultiEvalCudaDetail::HostDevPtr<OpMultiEvalCudaDetail::gridnode_t> node;
  OpMultiEvalCudaDetail::HostDevPtr<double> alpha;
  OpMultiEvalCudaDetail::HostDevPtr<double> pos;
  OpMultiEvalCudaDetail::HostDevPtr<double> data;
  OpMultiEvalCudaDetail::HostDevPtr<OpMultiEvalCudaDetail::limit_t> streamlimit;
  OpMultiEvalCudaDetail::HostDevPtr<uint32_t> levellimit;
  OpMultiEvalCudaDetail::HostDevPtr<uint32_t> subs;
  uint32_t maxlevel;

  uint32_t DIM;
  uint32_t N;
  uint32_t _N;
  uint32_t M;
  uint32_t _M;
  uint32_t polygrad;

 public:
  /// Constructor. Autocheck dataset for Morton order and turn on permutation if not
  OperationMultiEvalCuda(base::Grid& grid, base::DataMatrix& dataset, const uint32_t& grad);
  /// Constructor. Autocheck for Morton order can be disabled
  OperationMultiEvalCuda(base::Grid& grid, base::DataMatrix& dataset, const uint32_t& grad,
                         bool autocheck);
  /// Destructor
  ~OperationMultiEvalCuda();
  /// Standard evaluation
  void mult(base::DataVector& source, base::DataVector& result) override;
  /// Transposed evaluation
  void multTranspose(base::DataVector& source, base::DataVector& result) override;

  /// Does all preprocessing for given grid and dataset and copies data to the GPU
  void prepare() override;

  /// Returns time in s of last mult, multTransposed, multTransposedFMA
  double getDuration() override;

  /// Transposed evaluation with additional FMA. result = 1/M * B*source + lambda * prev
  void multTransposeFMA(base::DataVector& source, base::DataVector& prev, double lambda,
                        base::DataVector& result);
};

}  // namespace datadriven
}  // namespace sgpp

#endif  // OPERATIONMULTIEVALCUDA_HPP
