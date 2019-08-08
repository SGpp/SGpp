// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>

#include <chrono>
#include <vector>
#include <string>

#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/OperationPruneGraphOCL.hpp>
#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/KernelPruneGraph.hpp>

namespace sgpp {
namespace datadriven {
namespace DensityOCLMultiPlatform {

/// Operation for density based graph pruning
template<typename T>
class OperationPruneGraphOCLMultiPlatform : public OperationPruneGraphOCL {
 private:
  size_t dims;
  size_t gridSize;
  bool verbose;
  size_t dataSize;
  /// OpenCL kernel which executes the graph creation
  sgpp::datadriven::DensityOCLMultiPlatform::KernelPruneGraph<T> *graph_kernel;
  /// Vector with all OpenCL devices
  std::vector<std::shared_ptr<base::OCLDevice>> devices;
  /// OpenCL Manager
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;

  /// Copy of the grid used for the density estimation
  std::vector<int> pointsVector;
  /// Copy of the alpha vector used for the density estimation
  std::vector<T> alphaVector;
  /// Copy of the dataset
  std::vector<T> dataVector;

  size_t start_id;
  size_t chunksize;

 public:
  /// Constructor using a DataMatrix
  OperationPruneGraphOCLMultiPlatform(base::Grid& grid, base::DataVector& alpha,
                                      base::DataMatrix &data, size_t dims,
                                      std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                                      sgpp::base::OCLOperationConfiguration *parameters,
                                      T treshold, size_t k, size_t platform_id,
                                      size_t device_id) :
      OperationPruneGraphOCL(), dims(dims), gridSize(grid.getStorage().getSize()), verbose(false),
      dataSize(data.getSize()),  devices(manager->getDevices()),
      manager(manager), alphaVector(gridSize), dataVector(data.getSize()) {
    // Store Grid in a opencl compatible buffer
    sgpp::base::GridStorage& gridStorage = grid.getStorage();
    size_t pointscount = 0;
    for (size_t i = 0; i < gridSize; i++) {
      sgpp::base::HashGridPoint &point = gridStorage.getPoint(i);
      pointscount++;
      for (size_t d = 0; d < dims; d++) {
        pointsVector.push_back(point.getIndex(d));
        pointsVector.push_back(point.getLevel(d));
      }
    }
    for (size_t i = 0; i < gridSize; i++)
      alphaVector[i] = static_cast<T>(alpha[i]);
    double *data_raw = data.getPointer();
    for (size_t i = 0; i < data.getSize(); i++)
      dataVector[i] = static_cast<T>(data_raw[i]);
    // look for the chosen platform and device and create kernel with it
    size_t platformcounter = 0;
    size_t devicecounter = 0;
    cl_device_id old_device_id = devices[0]->deviceId;
    cl_platform_id old_platform_id = devices[0]->platformId;
    size_t counter = 0;
    bool success = false;
    for (auto device : devices) {
      if (device->platformId != old_platform_id) {
        platformcounter++;
        old_platform_id = device->platformId;
        devicecounter = 0;
      }
      if (device->deviceId != old_device_id) {
        devicecounter++;
        old_device_id = device->deviceId;
      }
      if (platformcounter == platform_id &&
          devicecounter == device_id) {
        json::Node &deviceNode =
            (*parameters)["PLATFORMS"][device->platformName]["DEVICES"][device->deviceName];
        json::Node &configuration = deviceNode["KERNELS"]["removeEdges"];
        graph_kernel = new KernelPruneGraph<T>(devices[counter], dims, treshold, k, manager,
                                               configuration, pointsVector, alphaVector,
                                               dataVector);
        if (configuration["VERBOSE"].getBool())
          verbose = true;
        success = true;
        break;
      }
      counter++;
    }

    // Check whether a kernel was created or not
    if (!success) {
      std::stringstream errorString;
      errorString << "OCL Error: Platform with index " << platform_id
                  << " and the device with index " << device_id << std::endl
                  << " not found! Please check your OpenCL installation!" << std::endl;
      throw base::operation_exception(errorString.str());
    }
  }
  /// Constructor using a double vector as a dataset and a serialized grid
  OperationPruneGraphOCLMultiPlatform(int *gridpoints, size_t gridSize, size_t dimensions,
                                      double *alpha, base::DataMatrix &data,
                                      std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                                      sgpp::base::OCLOperationConfiguration *parameters,
                                      T treshold, size_t k, size_t platform_id,
                                      size_t device_id) :
      OperationPruneGraphOCL(), dims(dimensions), gridSize(gridSize), verbose(false),
      dataSize(data.getSize()),  devices(manager->getDevices()),
      manager(manager), alphaVector(gridSize), dataVector(data.getSize()) {
    // Store Grid in a opencl compatible buffer
    std::vector<int> points;
    for (size_t i = 0; i < gridSize; i++) {
      for (size_t d = 0; d < dims; d++) {
        points.push_back(gridpoints[2 * dimensions * i + 2 * d]);
        points.push_back(gridpoints[2 * dimensions * i + 2 * d + 1]);
      }
    }
    if (verbose)
      std::cout << "Grid stored into integer array! Number of gridpoints: "
                << gridSize << std::endl;
    for (size_t i = 0; i < gridSize; i++)
      alphaVector[i] = static_cast<T>(alpha[i]);
    double *data_raw = data.getPointer();
    for (size_t i = 0; i < data.getSize(); i++)
      dataVector[i] = static_cast<T>(data_raw[i]);
    if (verbose)
      std::cout << "Data stored into float array! Number of datapoints: "
                << data.getSize() << std::endl;

    // look for the chosen platform and device and create kernel with it
    size_t platformcounter = 0;
    size_t devicecounter = 0;
    cl_device_id old_device_id = devices[0]->deviceId;
    cl_platform_id old_platform_id = devices[0]->platformId;
    size_t counter = 0;
    bool success = false;
    for (auto device : devices) {
      if (device->platformId != old_platform_id) {
        platformcounter++;
        old_platform_id = device->platformId;
        devicecounter = 0;
      }
      if (device->deviceId != old_device_id) {
        devicecounter++;
        old_device_id = device->deviceId;
      }
      if (platformcounter == platform_id &&
          devicecounter == device_id) {
        json::Node &deviceNode =
            (*parameters)["PLATFORMS"][device->platformName]["DEVICES"][device->deviceName];
        json::Node &configuration = deviceNode["KERNELS"]["removeEdges"];
        graph_kernel = new KernelPruneGraph<T>(devices[counter], dims, treshold, k, manager,
                                               configuration, points, alphaVector,
                                               dataVector);
        if (configuration["VERBOSE"].getBool())
          verbose = true;
        success = true;
        break;
      }
      counter++;
    }

    // Check whether a kernel was created or not
    if (!success) {
      std::stringstream errorString;
      errorString << "OCL Error: Platform with index " << platform_id
                  << " and the device with index " << device_id << std::endl
                  << " not found! Please check your OpenCL installation!" << std::endl;
      throw base::operation_exception(errorString.str());
    }
  }

  ~OperationPruneGraphOCLMultiPlatform() {
    delete graph_kernel;
  }

  /// Deletes all nodes and edges within areas of low density which are in the given graph chunk
  virtual void prune_graph(std::vector<int> &graph, size_t startid = 0, size_t chunksize = 0) {
    if (verbose)
      std::cout << "Pruning graph for" << graph.size() << " nodes" << std::endl;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    this->graph_kernel->prune_graph(graph, startid, chunksize);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    if (verbose)
      std::cout << "duration prune graph: " << elapsed_seconds.count() << std::endl;
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
