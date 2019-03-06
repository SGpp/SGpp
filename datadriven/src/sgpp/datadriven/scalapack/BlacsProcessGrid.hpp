/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * BlacsProcessGrid.hpp
 *
 * Created on: Jan 14, 2019
 *     Author: Jan Schopohl
 */
#ifdef USE_SCALAPACK

#pragma once

#include <sgpp/datadriven/scalapack/blacs.hpp>

#include <atomic>
#include <cstddef>

namespace sgpp {
namespace datadriven {

/**
 * This class represents a BLACS process grid for use with ScaLAPACK.
 */
class BlacsProcessGrid {
 public:
  /**
   * Creates a square process grid of maximum size.
   * Given p processes, creates a sqrt(p) * sqrt(p) process grid.
   * Always call this method from *all* processes, as the init
   * method of a BLACS grid has to be called from all processes, otherwise a deadlock will occur.
   */
  BlacsProcessGrid();

  /**
   * Creates a BLACS process grid with a certain number of rows and columns.
   * There must be at least rows * columns processes available.
   * Always call this method from *all* processes, as the init
   * method of a BLACS grid has to be called from all processes, otherwise a deadlock will occur.
   *
   * @param rows
   * @param columns
   */
  BlacsProcessGrid(int rows, int columns);

  /**
   * Cannot be copied, otherwise errors with multiple calls to blacs_gridexit are possible
   */
  BlacsProcessGrid(const BlacsProcessGrid&) = delete;

  ~BlacsProcessGrid();

  /**
   * @returns the context handle of the BLACS context
   */
  int getContextHandle() const;

  /**
   * @returns Total number of rows of the grid
   */
  int getTotalRows() const;

  /**
   * @returns Total number of columns of the grid
   */
  int getTotalColumns() const;

  /**
   * @returns Row of the current process
   */
  int getCurrentRow() const;

  /**
   * @returns Column of the current process
   */
  int getCurrentColumn() const;

  /**
   * @returns Number of the current process
   */
  int getCurrentProcess() const;

  /**
   * @returns True if the current process is part of this grid, else false
   */
  bool isProcessInGrid() const;

  /**
   * Can only be called after BLACS initialization.
   * @returns the number of available mpi processes
   */
  static size_t availableProcesses();

  /**
   * Initialize BLACS, should only be called once.
   */
  static void initializeBlacs();

  /**
   * Exit BLACS, should only be called once.
   */
  static void exitBlacs();

 private:
  // flag to check if BLACS was initialized
  static bool blacsInitialized;

  // system context for use in gridinit
  static int systemContext;

  // total number of available processes
  static int numberOfProcesses;

  // BLACS context handle
  int ictxt;

  // process grid rows
  int rows;

  // process grid columns
  int columns;

  // describes if the current process is part of this grid
  bool partOfGrid;

  // current process number
  int mypnum;

  // row of the current process
  int myrow;

  // column of the current process
  int mycolumn;
};
}  // namespace datadriven
}  // namespace sgpp

#endif  // USE_SCALAPACK