
#ifdef USE_SCALAPACK

#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>
#include <sgpp/datadriven/scalapack/blacs.hpp>
#include <sgpp/datadriven/scalapack/scalapack.hpp>

#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <memory>
#include <vector>

using sgpp::datadriven::BlacsProcessGrid;
using sgpp::datadriven::DataMatrixDistributed;

void test() {
  std::shared_ptr<BlacsProcessGrid> grid = std::make_shared<BlacsProcessGrid>(1, 1);
  std::cout << "Context handle: " << grid->getContextHandle() << std::endl;
  std::cout << "grid: " << grid->getTotalRows() << " x " << grid->getTotalColumns() << std::endl;
  std::cout << "process: " << grid->getCurrentRow() << ", " << grid->getCurrentColumn()
            << std::endl;

  DataMatrixDistributed matrix(grid, 10, 10, 128, 64, 0.0);
  std::cout << "matrix initialized" << std::endl;

  std::cout << "local size: " << matrix.getLocalRows() << " x " << matrix.getLocalColumns()
            << std::endl;

  for (int i = 0; i < matrix.getLocalRows(); i++) {
    for (int j = 0; j < matrix.getLocalColumns(); j++) {
      std::cout << matrix(i, j);
    }
    std::cout << std::endl;
  }

  std::vector<double> test{0.0, 1.0, 3.0, 0.0};

  DataMatrixDistributed matrix2(test.data(), grid, 2, 2, 2, 1);

  int m = 2;
  int n = 2;
  int ia = 1;
  int ja = 1;
  int irprnt = 0;
  int icprnt = 0;
  std::string identifier = "A";
  int nout = 0;
  auto work = std::vector<double>(m * n * 4);
  // sgpp::datadriven::pdlaprnt_(m, n, matrix2.data(), ia, ja, matrix2.getDescriptor(), irprnt,
  // icprnt,
  //                            identifier.c_str(), nout, work.data());

  for (int i = 0; i < matrix2.getLocalRows(); i++) {
    for (int j = 0; j < matrix2.getLocalColumns(); j++) {
      std::cout << matrix2.getLocalPointer()[(i * 2) + j];
    }
    std::cout << std::endl;
  }

  sgpp::base::DataMatrix local = matrix2.toLocalDataMatrix();

  std::cout << "process " << grid->getCurrentProcess()
            << " converted local size: " << local.getNrows() << " x " << local.getNcols()
            << std::endl;

  for (int i = 0; i < matrix2.getLocalRows(); i++) {
    for (int j = 0; j < matrix2.getLocalColumns(); j++) {
      std::cout << matrix2.getLocalPointer()[(i * 2) + j];
    }
    std::cout << std::endl;
  }

  std::cout << local.toString() << std::endl;
}

#endif /* USE_SCALAPACK */

int main(int argc, char** argv) {
#ifdef USE_SCALAPACK

  BlacsProcessGrid::initializeBlacs();
  test();
  BlacsProcessGrid::exitBlacs();

  // generate matrices
  /*int m = 8;
  int n = 8;
  int k = 8;

  std::vector<double> A = std::vector<double>(m * k);
  std::vector<double> B = std::vector<double>(k * n);
  std::vector<double> C = std::vector<double>(m * n);

  int ictxt = 100;
  int rows = 1;
  int cols = 1;

  sl_init_(ictxt, rows, cols);

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int n_rows, n_cols, myrow, mycol;

  blacs_gridinfo_(ictxt, n_rows, n_cols, myrow, mycol);

  std::cout << "This process is in row " << myrow << " out of " << n_rows << " and col " << mycol
            << " out of " << n_cols << std::endl;

  char uplo = 'A';
  int ia = 1;
  int ja = 1;
  double alpha = 1.0;
  double beta = 0.0;
  int desca[DLEN_];
  desca[DTYPE_] = 1;
  desca[CTXT_] = ictxt;
  desca[M_] = m;
  desca[N_] = n;
  desca[MB_] = 8;
  desca[NB_] = 8;
  desca[RSRC_] = 0;
  desca[CSRC_] = 0;
  desca[LLD_] = k;

  pdlaset_(uplo, m, n, alpha, beta, A.data(), ia, ja, desca);

  MPI_Barrier(MPI_COMM_WORLD);

  auto work = std::vector<double>(m * n * 4);
  std::string identifier = "A";
  int irprnt = myrow;
  int icprnt = mycol;
  int nout = 0;

  // pdlaprnt_(m, n, A.data(), ia, ja, desca, irprnt, icprnt, identifier.c_str(), nout,
  // work.data())

  std::cout << "A: " << std::endl;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      std::cout << A[(i * n) + j];
    }
    std::cout << std::endl;
  }

  // Finalize the MPI environment.
  blacs_gridexit_(ictxt);

  int cont = 0;
  blacs_exit_(cont);*/

#endif /*USE_SCALAPACK*/
}
