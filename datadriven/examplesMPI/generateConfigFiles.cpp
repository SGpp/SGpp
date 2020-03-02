// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_MPI

#include <sgpp/base/tools/OperationConfiguration.hpp>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OperationCreateGraphOCLSingleDevice.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/MPIEnviroment.hpp>
#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/OpFactory.hpp>
#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

#include <iostream>
#include <sstream>
#include <string>

int main(int argc, char **argv) {
  if (argc == 4) {
    int packagesize_master = std::stoi(argv[1]);
    int number_compute_nodes = std::stoi(argv[2]);
    sgpp::base::OperationConfiguration conf =
        sgpp::datadriven::clusteringmpi::MPIEnviroment::createMPIConfiguration(
            packagesize_master, number_compute_nodes);
    conf.serialize(argv[3]);
    return 0;
  }
  if (argc == 6) {
    int packagesize_master = std::stoi(argv[1]);
    int number_leutnant_nodes = std::stoi(argv[2]);
    int packagesize_leutnant = std::stoi(argv[3]);
    int number_slave_nodes = std::stoi(argv[4]);
    sgpp::base::OperationConfiguration conf =
        sgpp::datadriven::clusteringmpi::MPIEnviroment::createMPIConfiguration(
            packagesize_master, number_leutnant_nodes, packagesize_leutnant,
            number_slave_nodes);
    conf.serialize(argv[5]);
    return 0;
  } else {
    std::cout << "Usage:" << std::endl;
    std::cout
        << "./generateConfigFiles <Packagesize> <Number of compute nodes> "
           "<Filename>"
        << std::endl
        << "OR " << std::endl
        << "./generateConfigFiles <Packagesize> <Number of leutnant nodes> "
        << "<Packagesize leutnant nodes> <Number of compute nodes per "
           "leutnant> <Filename>"
        << std::endl;
    return 0;
  }
  return 0;
}

#else
#include <iostream>
int main(int argc, char** argv) {
  std::cout << "error: build with MPI to enable this example" << std::endl;
  return 0;
}
#endif
