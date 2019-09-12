// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <unistd.h>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/MPIEnviroment.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationCreateGraphMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationCreatePrunedGraphMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationDensityMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationRhsMPI.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <sgpp/datadriven/operation/hash/OperationMPI/OperationGridBaseMPI.hpp>

#include <iostream>
#include <string>
#include <vector>

#include "sgpp/datadriven/tools/ARFFTools.hpp"

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
    std::cout << "Wrong arguments - consult source file for more informations!" << std::endl;
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
  std::chrono::time_point<std::chrono::high_resolution_clock> grid_creation_start,
      grid_creation_end;
  grid_creation_start = std::chrono::system_clock::now();

  sgpp::base::Grid *grid = sgpp::base::Grid::createLinearGrid(dim);
  sgpp::base::GridGenerator &gridGen = grid->getGenerator();
  std::cout << "level:" << level << std::endl;
  gridGen.regular(level);
  size_t gridsize = grid->getStorage().getSize();
  std::cerr << "Grid created! Number of grid points:     " << gridsize << std::endl;

  grid_creation_end = std::chrono::system_clock::now();
  if (rank == 0) {
    std::cout << "grid creation duration: "
              << std::chrono::duration_cast<std::chrono::seconds>(grid_creation_end -
                                                                  grid_creation_start)
                     .count()
              << "s" << std::endl;
  }

  sgpp::datadriven::clusteringmpi::OperationDensityMultMPI mult_op(*grid, 0.001, "MyOCLConf.cfg");
  sgpp::base::DataVector alpha(gridsize);
  sgpp::base::DataVector result(gridsize);
  alpha.setAll(1.0);

  std::chrono::time_point<std::chrono::high_resolution_clock> rhs_start, rhs_end;
  rhs_start = std::chrono::system_clock::now();

  // Create right hand side vector
  sgpp::base::DataVector rhs(gridsize);
  sgpp::datadriven::clusteringmpi::OperationDensityRhsMPI rhs_op(*grid, dataset, "MyOCLConf.cfg");
  rhs_op.generate_b(rhs);

  rhs_end = std::chrono::system_clock::now();
  if (rank == 0) {
    std::cout << "rhs creation duration: "
              << std::chrono::duration_cast<std::chrono::seconds>(rhs_end - rhs_start).count()
              << "s" << std::endl;
  }

  // Solve for alpha vector via CG solver
  std::chrono::time_point<std::chrono::high_resolution_clock> solver_start, solver_end;
  solver_start = std::chrono::system_clock::now();

  alpha.setAll(1.0);
  sgpp::solver::ConjugateGradients solver(1000, 0.001);
  solver.solve(mult_op, alpha, rhs, false, true);

  solver_end = std::chrono::system_clock::now();
  if (rank == 0) {
    std::cout << "solver duration: "
              << std::chrono::duration_cast<std::chrono::seconds>(solver_end - solver_start).count()
              << "s" << std::endl;
  }

  double max = alpha.max();
  double min = alpha.min();
  for (size_t i = 0; i < gridsize; i++) alpha[i] = alpha[i] * 1.0 / (max - min);

  // Create and prune knn graph
  std::chrono::time_point<std::chrono::high_resolution_clock> create_knn_start, create_knn_end;
  create_knn_start = std::chrono::system_clock::now();
  std::cout << "Graph Creation/Pruning:" << std::endl;
  sgpp::datadriven::clusteringmpi::OperationPrunedGraphCreationMPI graph_op(
      *grid, alpha, dataset, 12, 0.2, "MyOCLConf.cfg");
  std::vector<int> knn_graph;
  graph_op.create_graph(knn_graph);

  create_knn_end = std::chrono::system_clock::now();
  if (rank == 0) {
    std::cout << "create knn operation duration: "
              << std::chrono::duration_cast<std::chrono::seconds>(create_knn_end - create_knn_start)
                     .count()
              << "s" << std::endl;
  }

  std::chrono::time_point<std::chrono::high_resolution_clock> find_clusters_start,
      find_clusters_end;
  find_clusters_start = std::chrono::system_clock::now();

  std::vector<size_t> cluster_assignments =
      sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCL::find_clusters(knn_graph,
                                                                                        12);

  find_clusters_end = std::chrono::system_clock::now();
  if (rank == 0) {
    std::cout << "find clusters duration: "
              << std::chrono::duration_cast<std::chrono::seconds>(find_clusters_end -
                                                                  find_clusters_start)
                     .count()
              << "s" << std::endl;
  }

  std::ofstream out("cluster_erg.txt");
  for (size_t datapoint : cluster_assignments) {
    out << datapoint << " ";
  }

  // Calc time
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  std::cout << "finished computation at " << std::ctime(&end_time)
            << "elapsed time: " << elapsed_seconds.count() << "s\n";
  // Cleanup MPI enviroment
  sgpp::datadriven::clusteringmpi::MPIEnviroment::release();
  return 0;
}
