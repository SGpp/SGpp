// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/base/opencl/QueueLoadBalancer.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/OpFactory.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <chrono>
#include <vector>
#include <string>
namespace sgpp {
namespace datadriven {
namespace ClusteringOCL {

template<typename T>
class OperationClusteringOCL {
 private:
  bool verbose;
  std::string opencl_configuration;

 public:
  OperationClusteringOCL(std::string opencl_configuration, bool verbose)
      : verbose(verbose), opencl_configuration(opencl_configuration) {
  }

  std::vector<size_t> calculate_clusters(base::Grid *grid, sgpp::base::DataMatrix &dataset,
                                         double lambda, size_t k, double treshold) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    size_t dimension = dataset.getNcols();
    size_t gridsize = grid->getStorage()->size();
    base::DataVector alpha(gridsize);
    base::DataVector b(gridsize);
    DensityOCLMultiPlatform::OperationDensity *mult_operation;
    DensityOCLMultiPlatform::OperationCreateGraphOCL *graph_operation;
    DensityOCLMultiPlatform::OperationPruneGraphOCL *prune_operation;

    start = std::chrono::system_clock::now();
    sg::solver::ConjugateGradients *solver = new sg::solver::ConjugateGradients(1000, 0.001);
    mult_operation = createDensityOCLMultiPlatformConfigured(*grid, dimension, lambda,
                                                             opencl_configuration);
    if (verbose)
      std::cout << "Creating rhs..." << std::endl;
    mult_operation->generateb(dataset, b);

    if (verbose)
      std::cout << "Creating alpha..." << std::endl;
    solver->solve(*mult_operation, alpha, b, false, verbose);
    double max = alpha.max();
    double min = alpha.min();
    for (size_t i = 0; i < gridsize; i++)
      alpha[i] = alpha[i]*1.0/(max-min);

    if (verbose)
      std::cout << "Starting graph creation..." << std::endl;
    graph_operation = createNearestNeighborGraphConfigured(dataset, k, dimension,
                                                           opencl_configuration);
    std::vector<int> graph(dataset.getNrows()*k);
    graph_operation->create_graph(graph);

    if (verbose)
      std::cout << "Starting graph pruning..." << std::endl;
    prune_operation = pruneNearestNeighborGraphConfigured(*grid, dimension, alpha, dataset,
                                                          treshold, k, opencl_configuration);
    prune_operation->prune_graph(graph);

    std::vector<size_t> ret = DensityOCLMultiPlatform::
        OperationCreateGraphOCL::find_clusters(graph, k);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    if (verbose) {
      std::cout << "Time required for clustering: " << elapsed_seconds.count() << std::endl;
    }

    return ret;
  }
};
}  // namespace ClusteringOCL
}  // namespace datadriven
}  // namespace sgpp
