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

#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>

#ifdef USE_SCALAPACK
#include <mpi.h>
#endif /* USE_SCALAPACK */
#include <unistd.h>
#include <cmath>
#include <iostream>
#include <string>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/datadriven/scalapack/blacs.hpp>
#include <sgpp/datadriven/scalapack/scalapack.hpp>

namespace sgpp {
namespace datadriven {

int BlacsProcessGrid::systemContext = 0;
int BlacsProcessGrid::mypnum = 0;
int BlacsProcessGrid::numberOfProcesses = 0;
bool BlacsProcessGrid::blacsInitialized = false;

BlacsProcessGrid::BlacsProcessGrid(int rows, int columns) : rows(rows), columns(columns) {
#ifdef USE_SCALAPACK
  if (!blacsInitialized) {
    throw sgpp::base::application_exception("BLACS not initialized!");
  }

  if (rows <= 0 || columns <= 0) {
    this->rows = static_cast<int>(std::sqrt(numberOfProcesses));
    this->columns = static_cast<int>(std::sqrt(numberOfProcesses));
  }

  if (this->rows * this->columns > availableProcesses()) {
    throw sgpp::base::application_exception(
        "Not enough processes available to form BLACS process grid!");
  }

  int ignore = -1;
  int systemContext = -1;
  Cblacs_get(ignore, 0, systemContext);
  ictxt = systemContext;
  Cblacs_gridinit(ictxt, "R", this->rows, this->columns);

  Cblacs_gridinfo(ictxt, this->rows, this->columns, myrow, mycolumn);

  if (myrow >= 0 && mycolumn >= 0) {
    Cblacs_pnum(ictxt, myrow, mycolumn);
    this->partOfGrid = true;
  } else {
    this->partOfGrid = false;
  }
#else
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

BlacsProcessGrid::~BlacsProcessGrid() {
#ifdef USE_SCALAPACK
  if (partOfGrid) {
    // only exit the grid if this process is actually part of it
    Cblacs_gridexit(ictxt);
  }
#endif /* USE_SCALAPACK */
}

int BlacsProcessGrid::getContextHandle() const { return this->ictxt; }

int BlacsProcessGrid::getTotalRows() const { return this->rows; }

int BlacsProcessGrid::getTotalColumns() const { return this->columns; }

int BlacsProcessGrid::getCurrentRow() const { return this->myrow; }

int BlacsProcessGrid::getCurrentColumn() const { return this->mycolumn; }

int BlacsProcessGrid::getRowColumnIndex() const { return (myrow * columns) + mycolumn; }

int BlacsProcessGrid::getProcessesInGrid() const { return rows * columns; }

int BlacsProcessGrid::getCurrentProcess() { return mypnum; }

bool BlacsProcessGrid::isProcessInGrid() const { return this->partOfGrid; }

int BlacsProcessGrid::availableProcesses() {
  if (!blacsInitialized) {
    throw sgpp::base::application_exception("BLACS not initialized!");
  }
  return numberOfProcesses;
}

void BlacsProcessGrid::initializeBlacs() {
#ifdef USE_SCALAPACK
  // init BLACS and the MPI environment
  MPI_Init(nullptr, nullptr);
  Cblacs_pinfo(mypnum, numberOfProcesses);
  blacsInitialized = true;

  // get the name of the node
  char nodeName[MPI_MAX_PROCESSOR_NAME];
  int nameLength;
  MPI_Get_processor_name(nodeName, &nameLength);
  std::string node(nodeName, nameLength);

  std::cout << "Init BLACS and MPI on node " << node << std::endl;
#else
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

void BlacsProcessGrid::exitBlacs() {
#ifdef USE_SCALAPACK
  // Finalizes ScaLAPACK, BLACS and the MPI environment
  std::cout << "Exit BLACS and MPI" << std::endl;
  Cblacs_exit(1);
  MPI_Finalize();
#else
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

}  // namespace datadriven
}  // namespace sgpp