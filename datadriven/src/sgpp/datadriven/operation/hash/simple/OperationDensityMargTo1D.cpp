// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalize.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

void OperationDensityMargTo1D::computeMarginalizationIndices(std::vector<size_t>& dim_x,
                                                             size_t numDims,
                                                             std::vector<size_t>& margDims) {
  if (dim_x.size() == numDims) {
    margDims.clear();
    return;
  }

  size_t count = 0;
  for (size_t idim = 0; idim < numDims; idim++) {
    size_t i = 0;
    bool marginalize = true;
    while (i < dim_x.size() && marginalize) {
      marginalize &= dim_x[i] != idim;
      i++;
    }

    if (marginalize) {
      margDims.push_back(idim - count);
      count++;
    }
  }
}

void OperationDensityMargTo1D::margToDimX(base::DataVector* alpha, base::Grid*& grid_x,
                                          base::DataVector*& alpha_x, size_t dim_x) {
  size_t numDims = this->grid->getDimension();
  //  size_t count = 0;

  if ((numDims > 1) && (dim_x <= numDims - 1)) {
    std::vector<size_t> margDims = {dim_x};
    margToDimXs(alpha, grid_x, alpha_x, margDims);
    //    marg_next_dim(this->grid, alpha, grid_x, alpha_x, dims, dim_x, count);
  } else if (numDims <= 1) {
    throw base::operation_exception(
        "Error: grid dimension is not greater than one. Operation aborted!");
  } else {
    throw base::operation_exception("Error: dimension out of range. Operation aborted!");
  }
  return;
}

void OperationDensityMargTo1D::margToDimXs(base::DataVector* alpha, base::Grid*& grid_x,
                                           base::DataVector*& alpha_x, std::vector<size_t>& dim_x) {
  size_t numDims = this->grid->getDimension();

  if (numDims == dim_x.size()) {
    // no marginalization required, just copy the grid and the coefficient vector
    grid_x = grid->clone();
    alpha_x->resize(alpha->getSize());
    for (size_t i = 0; i < alpha->getSize(); i++) {
      alpha_x->set(i, alpha->get(i));
    }
    return;
  }

  // prepare dimensions over which we want to integrate
  std::vector<size_t> margDims;
  computeMarginalizationIndices(dim_x, numDims, margDims);
  // do the integration
  marg_next_dim(grid, alpha, grid_x, alpha_x, margDims, 0);
}

void OperationDensityMargTo1D::marg_next_dim(base::Grid* g_in, base::DataVector* a_in,
                                             base::Grid*& g_out, base::DataVector*& a_out,
                                             std::vector<size_t> margDims, size_t ix) {
  unsigned int op_dim = static_cast<unsigned int>(margDims[ix]);
  base::Grid* g_tmp = nullptr;
  base::DataVector* a_tmp = new base::DataVector(1);
  op_factory::createOperationDensityMarginalize(*g_in)->doMarginalize(*a_in, g_tmp, *a_tmp, op_dim);

  if (ix + 1 < margDims.size()) {
    marg_next_dim(g_tmp, a_tmp, g_out, a_out, margDims, ix + 1);
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
