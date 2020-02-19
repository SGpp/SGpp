// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OpFactory.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/globaldef.hpp>
#include <string>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/KernelCreateGraph.hpp>
namespace sgpp {
namespace datadriven {

DensityOCLMultiPlatform::OperationCreateGraphOCL *createNearestNeighborGraphConfigured(
    base::DataMatrix &dataset, size_t k, size_t dimensions, std::string opencl_conf,
    size_t platformid, size_t deviceid) {
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;

  std::cout << "Using configuration file " << opencl_conf << std::endl;
  sgpp::base::OCLOperationConfiguration *parameters =
      new sgpp::base::OCLOperationConfiguration(opencl_conf);
  manager = std::make_shared<base::OCLManagerMultiPlatform>((*parameters)["VERBOSE"].getBool());
  DensityOCLMultiPlatform::OperationCreateGraphOCL::load_default_parameters(parameters);

  if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
    return new DensityOCLMultiPlatform::OperationCreateGraphOCLSingleDevice<float>(
        dataset, dimensions, manager, parameters, k, platformid, deviceid);
  } else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
    return new DensityOCLMultiPlatform::OperationCreateGraphOCLSingleDevice<double>(
        dataset, dimensions, manager, parameters, k, platformid, deviceid);
  } else {
    std::stringstream errorString;
    errorString << "Error creating operation\"CreateGraphOCL\": "
                << " invalid value for parameter \"INTERNAL_PRECISION\"";
    throw base::factory_exception(errorString.str().c_str());
  }
  return nullptr;
}
DensityOCLMultiPlatform::OperationCreateGraphOCL *createNearestNeighborGraphConfigured(
    double *dataset, size_t dataset_size, size_t k, size_t dimensions,
    sgpp::base::OCLOperationConfiguration *parameters, size_t platformid, size_t deviceid) {
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;
  DensityOCLMultiPlatform::OperationCreateGraphOCL::load_default_parameters(parameters);
  manager = std::make_shared<base::OCLManagerMultiPlatform>((*parameters)["VERBOSE"].getBool());

  if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
    return new DensityOCLMultiPlatform::OperationCreateGraphOCLSingleDevice<float>(
        dataset, dataset_size, dimensions, manager, parameters, k, platformid, deviceid);
  } else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
    return new DensityOCLMultiPlatform::OperationCreateGraphOCLSingleDevice<double>(
        dataset, dataset_size, dimensions, manager, parameters, k, platformid, deviceid);
  } else {
    std::stringstream errorString;
    errorString << "Error creating operation\"CreateGraphOCL\": "
                << " invalid value for parameter \"INTERNAL_PRECISION\"";
    throw base::factory_exception(errorString.str().c_str());
  }
  return nullptr;
}

DensityOCLMultiPlatform::OperationCreateGraphOCL *createNearestNeighborGraphConfigured(
    double *dataset, size_t dataset_size, size_t k, size_t dimensions, std::string opencl_conf,
    size_t platformid, size_t deviceid) {
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;

  std::cout << "Using configuration file " << opencl_conf << std::endl;
  sgpp::base::OCLOperationConfiguration *parameters =
      new sgpp::base::OCLOperationConfiguration(opencl_conf);
  manager = std::make_shared<base::OCLManagerMultiPlatform>((*parameters)["VERBOSE"].getBool());
  DensityOCLMultiPlatform::OperationCreateGraphOCL::load_default_parameters(parameters);

  if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
    return new DensityOCLMultiPlatform::OperationCreateGraphOCLSingleDevice<float>(
        dataset, dataset_size, dimensions, manager, parameters, k, platformid, deviceid);
  } else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
    return new DensityOCLMultiPlatform::OperationCreateGraphOCLSingleDevice<double>(
        dataset, dataset_size, dimensions, manager, parameters, k, platformid, deviceid);
  } else {
    std::stringstream errorString;
    errorString << "Error creating operation\"CreateGraphOCL\": "
                << " invalid value for parameter \"INTERNAL_PRECISION\"";
    throw base::factory_exception(errorString.str().c_str());
  }
  return nullptr;
}
DensityOCLMultiPlatform::OperationCreateGraphOCL *createNearestNeighborGraphConfigured(
    base::DataMatrix &dataset, size_t k, size_t dimensions, std::string opencl_conf) {
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;

  std::cout << "Using configuration file " << opencl_conf << std::endl;
  sgpp::base::OCLOperationConfiguration *parameters =
      new sgpp::base::OCLOperationConfiguration(opencl_conf);
  manager = std::make_shared<base::OCLManagerMultiPlatform>((*parameters)["VERBOSE"].getBool());
  DensityOCLMultiPlatform::OperationCreateGraphOCL::load_default_parameters(parameters);

  size_t platformid = 0;
  if (parameters->contains("USE_PLATFORM") == true) {
    platformid = (*parameters)["USE_PLATFORM"].getInt();
  } else {
    std::stringstream errorString;
    errorString << "Error creating operation\"CreateGraphOCL\": "
                << "There is no given information about which opencl platform (ID) to use!"
                << " Add \"USE_PLATFORM\": platform_id to your configuration file "
                << "or use a different factory method." << std::endl;
    throw base::factory_exception(errorString.str().c_str());
  }
  size_t deviceid = 0;
  if (parameters->contains("USE_DEVICE") == true) {
    deviceid = (*parameters)["USE_DEVICE"].getInt();
  } else {
    std::stringstream errorString;
    errorString << "Error creating operation\"CreateGraphOCL\": "
                << "There is no given information about which opencl device (ID) to use!"
                << " Add \"USE_DEVICE\": device_id to your configuration file "
                << "or use a different factory method." << std::endl;
    throw base::factory_exception(errorString.str().c_str());
  }

  if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
    return new DensityOCLMultiPlatform::OperationCreateGraphOCLSingleDevice<float>(
        dataset, dimensions, manager, parameters, k, platformid, deviceid);
  } else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
    return new DensityOCLMultiPlatform::OperationCreateGraphOCLSingleDevice<double>(
        dataset, dimensions, manager, parameters, k, platformid, deviceid);
  } else {
    std::stringstream errorString;
    errorString << "Error creating operation\"CreateGraphOCL\": "
                << " invalid value for parameter \"INTERNAL_PRECISION\"";
    throw base::factory_exception(errorString.str().c_str());
  }
  return nullptr;
}
}  // namespace datadriven
}  // namespace sgpp
