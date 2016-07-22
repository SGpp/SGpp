// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/OperationEvalSGKernel/OperationEvalSGKernel.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace datadriven {

//Sparse Grid Kernel
// -------------------- constructors and destructors --------------------
OperationEvalSGKernel::OperationEvalSGKernel(sgpp::base::Grid& grid)
	: grid(grid) {}

OperationEvalSGKernel::~OperationEvalSGKernel() {}
// ----------------------------------------------------------------------


void OperationEvalSGKernel::phi(sgpp::base::DataVector& x, sgpp::base::DataVector& xTrans, size_t dataDim) {
  base::DataMatrix xMatrix(1,dataDim);
  xMatrix.setRow(0,x);
  base::DataVector alpha(1,1.0);
  op_factory::createOperationMultipleEval(grid, xMatrix)->multTranspose(alpha, xTrans);
}

}  // namespace datadriven
}  // namespace sgpp
