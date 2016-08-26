// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <sstream>
#include <string>

#include <sgpp/base/tools/OperationConfiguration.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/MPIEnviroment.hpp>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OperationCreateGraphOCLSingleDevice.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/OpFactory.hpp>
#include "sgpp/base/opencl/OCLManagerMultiPlatform.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"


int main(int argc, char **argv) {
  if (argc != 6) {
    std::cout << "Usage:" << std::endl;
    std::cout << "./generateConfigFiles <OCL config filename> <MPI config filename> "
              << "<Number of compute nodes> <Nodes packagesize> <devices packagesize>"
              << std::endl;
    return 0;
  }
  // Create opencl config file
  sgpp::base::OCLManagerMultiPlatform manager(true);
  auto configuration = manager.getConfiguration();
  sgpp::datadriven::DensityOCLMultiPlatform::
      OperationDensityOCL::load_default_parameters(configuration.get());
  sgpp::datadriven::DensityOCLMultiPlatform::
      OperationCreateGraphOCL::load_default_parameters(configuration.get());
  sgpp::datadriven::DensityOCLMultiPlatform::
      OperationPruneGraphOCL::load_default_parameters(configuration.get());
  configuration->serialize(argv[1]);

  // Get arguments
  size_t compute_nodes = 1;
  bool verbose = false, prefetching = true;
  int compute_node_packagesize = 20000;
  int device_packagesize = 4000;
  compute_nodes = std::stoi(argv[3]);
  compute_node_packagesize = std::stoi(argv[4]);
  device_packagesize = std::stoi(argv[5]);

  // Create MPI config file
  auto opencl_devices = manager.getDevices();
  sgpp::base::OperationConfiguration conf;
  conf.addIDAttr("VERBOSE", verbose);
  conf.addIDAttr("PACKAGE_SIZE", UINT64_C(10240));
  conf["PACKAGE_SIZE"].setInt(compute_node_packagesize);
  conf.addIDAttr("PREFETCHING", prefetching);
  std::unique_ptr<json::Node> workers(new json::DictNode);
  for (size_t i = 0; i < compute_nodes; ++i) {
    std::unique_ptr<json::Node> node_worker(new json::DictNode);
    node_worker->addIDAttr("VERBOSE", verbose);
    node_worker->addIDAttr("PACKAGE_SIZE", UINT64_C(2560));
    (*node_worker)["PACKAGE_SIZE"].setInt(device_packagesize);
    node_worker->addIDAttr("PREFETCHING", prefetching);
    std::unique_ptr<json::Node> device_workers(new json::DictNode);
    int node_counter = 0;
    int platform_counter = 0;
    int device_counter = 0;
    cl_device_id old_device_id = opencl_devices[0]->deviceId;
    cl_platform_id old_platform_id = opencl_devices[0]->platformId;
    for (std::shared_ptr<sgpp::base::OCLDevice> device : opencl_devices) {
      if (device->platformId != old_platform_id) {
        platform_counter++;
        old_platform_id = device->platformId;
      }
      if (device->deviceId != old_device_id) {
        device_counter++;
        old_device_id = device->deviceId;
      }
      std::cout << "device: " << device->platformId << " " << device->deviceId << std::endl;
      std::cout << "device: " << device->platformName << " " << device->deviceName << std::endl;
        std::unique_ptr<json::Node> device_worker(new json::DictNode);
        device_worker->addIDAttr("VERBOSE", verbose);
        device_worker->addIDAttr("OPENCL_PLATFORM", UINT64_C(0));
        (*device_worker)["OPENCL_PLATFORM"].setInt(platform_counter);
        device_worker->addIDAttr("OPENCL_DEVICE", UINT64_C(0));
        (*device_worker)["OPENCL_DEVICE"].setInt(device_counter);
        device_worker->addIDAttr("PREFETCHING", prefetching);
        std::string id = std::string("OPENCL_WORKER_") + std::to_string(node_counter);
        device_workers->addAttribute(id, std::move(device_worker));
        node_counter++;
    }
    node_worker->addAttribute("SLAVES", std::move(device_workers));
    std::string id = std::string("COMPUTE_NODE_") + std::to_string(i);
    workers->addAttribute(id, std::move(node_worker));
  }
  conf.addAttribute("SLAVES", std::move(workers));
  conf.serialize(argv[2]);
  return 0;
}
