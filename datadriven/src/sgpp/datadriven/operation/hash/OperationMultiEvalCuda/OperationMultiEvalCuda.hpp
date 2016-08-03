#ifndef OPERATIONMULTIEVALCUDA_HPP
#define OPERATIONMULTIEVALCUDA_HPP

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>

#include "cudaHelper.hpp"
#include "basicCuda.hpp"
#include <stdint.h>

namespace sgpp {
namespace datadriven {

/// OperationMultipleEval for polynomial basis functions (grad >= 2) using CUDA on grids without boundary nodes
class OperationMultiEvalCuda: public base::OperationMultipleEval {
 protected:
  sgpp::base::SGppStopwatch myTimer;
  double duration;
 
  HostDevPtr<gridnode_t> node;
  HostDevPtr<double> alpha;
  HostDevPtr<double> pos;
  HostDevPtr<double> data;
  HostDevPtr<limit_t> streamlimit;
  HostDevPtr<uint32_t> levellimit;
  HostDevPtr<uint32_t> subs;
  uint32_t maxlevel;
  
  uint32_t DIM;
  uint32_t N;
  uint32_t _N;
  uint32_t M;
  const bool _ordered;
  uint32_t polygrad;
 public:
  /// Constructor for datasets without a specific order
  OperationMultiEvalCuda(base::Grid& grid, base::DataMatrix& dataset, const uint32_t& grad);
  /// Constructor for datasets that are aligned along a Morton order curve
  OperationMultiEvalCuda(base::Grid& grid, base::DataMatrix& dataset, const uint32_t& grad, bool ordered);
  /// Destructor
  ~OperationMultiEvalCuda();
  /// Standard evaluation
  void mult(sgpp::base::DataVector& source, sgpp::base::DataVector& result) override;
  /// Transposed evaluation
  void multTranspose(sgpp::base::DataVector& source, sgpp::base::DataVector& result) override;
  
  /// Does all preprocessing for given grid and dataset and copies data to the GPU
  void prepare() override;
  
  /// Returns time in s of last mult, multTransposed, multTransposedFMA
  double getDuration() override;
  
  /// Transposed evaluation with additional FMA. result = 1/M * B*source + lambda * prev
  void multTransposeFMA(sgpp::base::DataVector& source, sgpp::base::DataVector& prev, double lambda, sgpp::base::DataVector& result);
};

}  // datadriven
}  // sgpp

#endif // OPERATIONMULTIEVALCUDA_HPP
