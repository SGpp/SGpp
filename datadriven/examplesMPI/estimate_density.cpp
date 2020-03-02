// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_MPI

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/MPIEnviroment.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationCreateGraphMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationCreatePrunedGraphMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationDensityMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationGridBaseMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationRhsMPI.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <unistd.h>

#include <iostream>
#include <string>
#include <vector>

// arg 1: MPI Config File
// arg 2: Dataset
// arg 3: Gridlevel
int main(int argc, char *argv[]) {
  // Init MPI enviroment - always has to be done first - capture slaves
  sgpp::datadriven::clusteringmpi::MPIEnviroment::init(argc, argv, true);

  // Measure times
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  start = std::chrono::system_clock::now();

  if (argc != 4) {
    std::cout << "Wrong arguments - consult source file for more informations!"
              << std::endl;
    sgpp::datadriven::clusteringmpi::MPIEnviroment::release();
    return 0;
  }

  // Loading dataset
  std::string filename = argv[2];
  std::cout << "Loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset data =
      sgpp::datadriven::ARFFTools::readARFFFromFile(filename);
  sgpp::base::DataMatrix &dataset = data.getData();
  size_t dim = data.getDimension();
  size_t level = std::stoi(argv[3]);

  sgpp::base::OperationConfiguration network_conf(argv[1]);
  sgpp::datadriven::clusteringmpi::MPIEnviroment::connect_nodes(network_conf);

  int rank = sgpp::datadriven::clusteringmpi::MPIEnviroment::get_node_rank();

  // Create Grid
  std::chrono::time_point<std::chrono::high_resolution_clock>
      grid_creation_start, grid_creation_end;
  grid_creation_start = std::chrono::system_clock::now();

  sgpp::base::Grid *grid = sgpp::base::Grid::createLinearGrid(dim);
  sgpp::base::GridGenerator &gridGen = grid->getGenerator();
  std::cout << "level:" << level << std::endl;
  gridGen.regular(level);
  size_t gridsize = grid->getStorage().getSize();
  std::cerr << "Grid created! Number of grid points:     " << gridsize
            << std::endl;

  grid_creation_end = std::chrono::system_clock::now();
  if (rank == 0) {
    std::cout << "grid creation duration: "
              << std::chrono::duration_cast<std::chrono::seconds>(
                     grid_creation_end - grid_creation_start)
                     .count()
              << "s" << std::endl;
  }

  sgpp::datadriven::clusteringmpi::OperationDensityMultMPI mult_op(
      *grid, 0.001, "MyOCLConf.cfg");
  sgpp::base::DataVector alpha(gridsize, 1.0);
  sgpp::base::DataVector result(gridsize);

  std::chrono::time_point<std::chrono::high_resolution_clock> rhs_start,
      rhs_end;
  rhs_start = std::chrono::system_clock::now();

  // Create right hand side vector
  sgpp::base::DataVector rhs(gridsize);
  sgpp::datadriven::clusteringmpi::OperationDensityRhsMPI rhs_op(
      *grid, dataset, "MyOCLConf.cfg");
  rhs_op.generate_b(rhs);

  rhs_end = std::chrono::system_clock::now();
  if (rank == 0) {
    std::cout << "rhs creation duration: "
              << std::chrono::duration_cast<std::chrono::seconds>(rhs_end -
                                                                  rhs_start)
                     .count()
              << "s" << std::endl;
  }

  // Solve for alpha vector via CG solver
  std::chrono::time_point<std::chrono::high_resolution_clock> solver_start,
      solver_end;
  solver_start = std::chrono::system_clock::now();

  alpha.setAll(1.0);
  sgpp::solver::ConjugateGradients solver(1000, 0.001);
  solver.solve(mult_op, alpha, rhs, false, true);

  solver_end = std::chrono::system_clock::now();
  if (rank == 0) {
    std::cout << "solver duration: "
              << std::chrono::duration_cast<std::chrono::seconds>(solver_end -
                                                                  solver_start)
                     .count()
              << "s" << std::endl;
  }

  // Cleanup MPI enviroment
  sgpp::datadriven::clusteringmpi::MPIEnviroment::release();
  return 0;
}

#else
#include <iostream>
int main(int argc, char** argv) {
  std::cout << "error: build with MPI to enable this example" << std::endl;
  return 0;
}
#endif
