// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1DLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalize.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {
void OperationDensityMargTo1DLinear::margToDimX(base::DataVector* alpha, base::Grid*& grid_x,
                                                base::DataVector*& alpha_x, size_t dim_x) {
  size_t dims = this->grid->getDimension();
  size_t count = 0;

  if ((dims > 1) && (dim_x <= dims - 1))
    marg_next_dim(this->grid, alpha, grid_x, alpha_x, dims, dim_x, count);
  else if (dims <= 1)
    throw base::operation_exception(
        "Error: grid dimension is not greater than one. Operation aborted!");
  else
    throw base::operation_exception("Error: dimension out of range. Operation aborted!");

  return;
}

void OperationDensityMargTo1DLinear::marg_next_dim(base::Grid* g_in, base::DataVector* a_in,
                                                   base::Grid*& g_out, base::DataVector*& a_out,
                                                   size_t dims, size_t dim_x, size_t& count) {
  unsigned int op_dim = (count < dim_x) ? 0 : 1;

  base::Grid* g_tmp = NULL;
  base::DataVector* a_tmp = new base::DataVector(1);
  op_factory::createOperationDensityMarginalize(*g_in)->doMarginalize(*a_in, g_tmp, *a_tmp, op_dim);

  count++;

  if (g_tmp->getDimension() > 1) {
    marg_next_dim(g_tmp, a_tmp, g_out, a_out, dims, dim_x, count);
    delete g_tmp;
    delete a_tmp;
  } else {
    g_out = g_tmp;
    a_out = a_tmp;
  }

  return;
}
}  // namespace datadriven
}  // namespace sgpp
