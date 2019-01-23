/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DistributedDataMatrix.cpp
 *
 * Created on: Jan 14, 2019
 *     Author: Jan Schopohl
 */
#include <sgpp/datadriven/scalapack/DistributedDataMatrix.hpp>

#include <sgpp/base/exception/not_implemented_exception.hpp>

namespace sgpp {
namespace datadriven {

DistributedDataMatrix::DistributedDataMatrix(std::shared_ptr<BlacsProcessGrid> grid, int globalRows,
                                             int globalColumns, DistributedDataMatrix::DTYPE dtype,
                                             int columnBlockSize, int rowBlockSize, double value)
    : DataMatrix(), grid(grid), globalRows(globalRows), globalColumns(globalColumns) {
  descriptor[dtype_] = static_cast<int>(dtype);

  // TODO(jan) implement other matrix types
  if (dtype != DistributedDataMatrix::DTYPE::DENSE) {
    throw sgpp::base::not_implemented_exception("Only Dense Matrix implemented at the moment");
  }

  localRows = numroc_(globalRows, rowBlockSize, grid->getCurrentRow(), 0, grid->getTotalRows());
  localColumns =
      numroc_(globalColumns, columnBlockSize, grid->getCurrentColumn(), 0, grid->getTotalColumns());

  // initialize local matrix
  this->assign(localRows * localColumns, value);

  descriptor[ctxt_] = grid->getContextHandle();
  descriptor[m_] = globalRows;
  descriptor[n_] = globalColumns;
  descriptor[mb_] = rowBlockSize;
  descriptor[nb_] = columnBlockSize;
  descriptor[rsrc_] = 0;
  descriptor[csrc_] = 0;
  descriptor[lld_] = localRows;  // TODO(jan) Can we use this to get row-major matrices?
}

DistributedDataMatrix::DistributedDataMatrix(double* input, std::shared_ptr<BlacsProcessGrid> grid,
                                             int globalRows, int globalColumns,
                                             DistributedDataMatrix::DTYPE dtype,
                                             int columnBlockSize, int rowBlockSize)
    : DistributedDataMatrix(grid, globalRows, globalColumns, dtype, columnBlockSize, rowBlockSize,
                            0.0) {
  // initialize local matrix A on process 0
  BlacsProcessGrid localGrid = BlacsProcessGrid(1, 1);
  int desca[dlen_];
  desca[ctxt_] = localGrid.getContextHandle();
  desca[m_] = globalRows;
  desca[n_] = globalColumns;
  desca[mb_] = globalRows;
  desca[nb_] = globalColumns;
  desca[rsrc_] = 0;
  desca[csrc_] = 0;
  descriptor[lld_] = globalRows;

  Cpdgemr2d(globalRows, globalColumns, input, 0, 0, desca, this->getPointer(), 0, 0,
            this->descriptor, this->grid->getContextHandle());
}

int* DistributedDataMatrix::getDescriptor() { return this->descriptor; }

int DistributedDataMatrix::getLocalRows() const { return this->localRows; }

int DistributedDataMatrix::getLocalColumns() const { return this->localColumns; }

}  // namespace datadriven
}  // namespace sgpp