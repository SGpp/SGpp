/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * BlacsProcessGrid.cpp
 *
 * Created on: Jan 14, 2019
 *     Author: Jan Schopohl
 */
#ifdef USE_SCALAPACK

#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>

#include <cmath>
#include <iostream>
#include <sgpp/datadriven/scalapack/blacs.hpp>
#include <sgpp/datadriven/scalapack/scalapack.hpp>

namespace sgpp {
namespace datadriven {

int BlacsProcessGrid::systemContext;
int BlacsProcessGrid::numberOfProcesses;

BlacsProcessGrid::BlacsProcessGrid() {
  size_t rows, columns = std::sqrt(numberOfProcesses);
  BlacsProcessGrid(rows, columns);
}

BlacsProcessGrid::BlacsProcessGrid(int rows, int columns) : rows(rows), columns(columns) {
  int ignore = -1;
  blacs_get_(ignore, 0, systemContext);
  ictxt = systemContext;
  blacs_gridinit_(ictxt, "R", this->rows, this->columns);

  // TODO(jan) error handling
  blacs_gridinfo_(ictxt, this->rows, this->columns, myrow, mycolumn);

  this->mypnum = blacs_pnum_(ictxt, myrow, mycolumn);
}

BlacsProcessGrid::~BlacsProcessGrid() {
  if (this->mypnum >= 0) {
    // only exit the grid if this process is actually part of it
    blacs_gridexit_(ictxt);
  }
}

int BlacsProcessGrid::getContextHandle() const { return this->ictxt; }

int BlacsProcessGrid::getTotalRows() const { return this->rows; }

int BlacsProcessGrid::getTotalColumns() const { return this->columns; }

int BlacsProcessGrid::getCurrentRow() const { return this->mycolumn; }

int BlacsProcessGrid::getCurrentColumn() const { return this->myrow; }

int BlacsProcessGrid::getCurrentProcess() const { return this->mypnum; }

void BlacsProcessGrid::initializeBlacs() {
  // init BLACS and the MPI environment
  int ignore = -1;
  blacs_get_(ignore, 0, systemContext);
  int mypnum;
  blacs_pinfo_(mypnum, numberOfProcesses);
}

void BlacsProcessGrid::exitBlacs() {
  // Finalizes ScaLAPACK, BLACS and the MPI environment
  blacs_exit_(0);
}

}  // namespace datadriven
}  // namespace sgpp

#endif  // USE_SCALAPACK