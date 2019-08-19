// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>

#include <chrono>
#include <vector>
#include <string>

#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OperationCreateGraphOCL.hpp>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/KernelCreateGraph.hpp>

namespace sgpp {
namespace datadriven {
namespace DensityOCLMultiPlatform {
/// Operation for creation of a k nearest neighbors graph using a dataset
template<typename T>
class OperationCreateGraphOCLSingleDevice : public OperationCreateGraphOCL {
 private:
  size_t dims;
  /// OpenCL kernel which executes the graph creation
  KernelCreateGraph<T> *graph_kernel;
  bool verbose;
  /// Vector with all OpenCL devices
  std::vector<std::shared_ptr<base::OCLDevice>> devices;
  /// OpenCL Manager
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;
  /// Copy of the dataset
  std::vector<T> dataVector;

 public:
  /// Constructor using a DataMatrix
  OperationCreateGraphOCLSingleDevice(base::DataMatrix& data, size_t dimensions,
                                      std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                                      sgpp::base::OCLOperationConfiguration *parameters,
                                      size_t k, size_t platform_id, size_t device_id) :
      OperationCreateGraphOCL(), dims(dimensions), verbose(false),
      devices(manager->getDevices()), manager(manager), dataVector(data.getSize()) {
    // put data into an vector with chosen precision
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
        json::Node &configuration = deviceNode["KERNELS"]["connectNeighbors"];
        graph_kernel = new KernelCreateGraph<T>(devices[counter], dims, k, dataVector,
                                                manager, configuration);
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
  /// Constructor which accepts double vector instead of a DataMatrix
  OperationCreateGraphOCLSingleDevice(double *dataset, size_t datasize, size_t dimensions,
                                      std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                                      sgpp::base::OCLOperationConfiguration *parameters,
                                      size_t k, size_t platform_id,
                                      size_t device_id) :
      OperationCreateGraphOCL(), dims(dimensions), verbose(false),
      devices(manager->getDevices()), manager(manager), dataVector(datasize) {
    // put data into an vector with chosen precision
    for (size_t i = 0; i < datasize; i++)
      dataVector[i] = static_cast<T>(dataset[i]);

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
        json::Node &configuration = deviceNode["KERNELS"]["connectNeighbors"];
        graph_kernel = new KernelCreateGraph<T>(devices[counter], dims, k, dataVector,
                                                manager, configuration);
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

  virtual ~OperationCreateGraphOCLSingleDevice(void) {
    delete graph_kernel;
  }
  void set_problemchunk(size_t start_id, size_t chunksize) {
    this->start_id = start_id;
    this->chunksize = chunksize;
  }

  /// Creates part of the k nearest neighbor graph (or all of it with the default parameters)
  void create_graph(std::vector<int> &resultVector, int startid = 0, int chunksize = 0) {
    if (verbose)
      std::cout << "Creating graph for " << dataVector.size() << " datapoints" << std::endl;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    try {
      graph_kernel->begin_graph_creation(startid, chunksize);
      graph_kernel->finalize_graph_creation(resultVector, startid, chunksize);
    }
    catch(base::operation_exception &e) {
      std::cerr << "Error! Could not create graph." << std::endl
                <<"Error Message: " << e.what() << std::endl;
      return;
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    if (verbose)
      std::cout << "duration create graph" << elapsed_seconds.count() << std::endl;
  }
  void begin_graph_creation(int startid, int chunksize) {
    graph_kernel->begin_graph_creation(startid, chunksize);
  }
  void finalize_graph_creation(std::vector<int> &resultVector, int startid, int chunksize) {
    graph_kernel->finalize_graph_creation(resultVector, startid, chunksize);
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
