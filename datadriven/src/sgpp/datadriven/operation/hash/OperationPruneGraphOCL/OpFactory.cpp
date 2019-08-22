// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/globaldef.hpp>
#include <string>
#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/OperationPruneGraphOCL.hpp>
namespace sgpp {
namespace datadriven {

DensityOCLMultiPlatform::OperationPruneGraphOCL*
pruneNearestNeighborGraphConfigured(base::Grid& grid, size_t dimensions, base::DataVector &alpha,
                                    base::DataMatrix &data, double treshold, size_t k,
                                    std::string opencl_conf, size_t platformid, size_t deviceid) {
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;

  std::cout << "Using configuration file " << opencl_conf << std::endl;
  sgpp::base::OCLOperationConfiguration *parameters =
      new sgpp::base::OCLOperationConfiguration(opencl_conf);
  manager = std::make_shared<base::OCLManagerMultiPlatform>((*parameters)["VERBOSE"].getBool());

  DensityOCLMultiPlatform::OperationPruneGraphOCL::load_default_parameters(parameters);

  if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
    return new DensityOCLMultiPlatform::
        OperationPruneGraphOCLMultiPlatform<float>(grid, alpha, data, dimensions, manager,
                                                   parameters,
                                                   static_cast<float>(treshold), k,
                                                   platformid, deviceid);
  } else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
    return new DensityOCLMultiPlatform::
        OperationPruneGraphOCLMultiPlatform<double>(grid, alpha, data, dimensions, manager,
                                                    parameters, treshold, k, platformid,
                                                    deviceid);
  } else {
    std::stringstream errorString;
    errorString << "Error creating operation\"OperationPruneGraphOCL\": "
                << " invalid value for parameter \"INTERNAL_PRECISION\"";
    throw base::factory_exception(errorString.str().c_str());
  }
  return nullptr;
}
DensityOCLMultiPlatform::OperationPruneGraphOCL*
pruneNearestNeighborGraphConfigured(int *gridpoints, size_t gridsize, size_t dimensions,
                                    double *alpha, base::DataMatrix &data,
                                    double treshold, size_t k, std::string opencl_conf,
                                    size_t platformid, size_t deviceid) {
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;

  std::cout << "Using configuration file " << opencl_conf << std::endl;
  sgpp::base::OCLOperationConfiguration *parameters =
      new sgpp::base::OCLOperationConfiguration(opencl_conf);
  manager = std::make_shared<base::OCLManagerMultiPlatform>((*parameters)["VERBOSE"].getBool());
  DensityOCLMultiPlatform::OperationPruneGraphOCL::load_default_parameters(parameters);

  if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
    return new DensityOCLMultiPlatform::
        OperationPruneGraphOCLMultiPlatform<float>(gridpoints, gridsize, dimensions, alpha, data,
                                                   manager, parameters,
                                                   static_cast<float>(treshold), k,
                                                   platformid, deviceid);
  } else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
    return new DensityOCLMultiPlatform::
        OperationPruneGraphOCLMultiPlatform<double>(gridpoints, gridsize, dimensions, alpha, data,
                                                    manager, parameters, treshold, k, platformid,
                                                    deviceid);
  } else {
    std::stringstream errorString;
    errorString << "Error creating operation\"OperationPruneGraphOCL\": "
                << " invalid value for parameter \"INTERNAL_PRECISION\"";
    throw base::factory_exception(errorString.str().c_str());
  }
  return nullptr;
}

DensityOCLMultiPlatform::OperationPruneGraphOCL*
pruneNearestNeighborGraphConfigured(int *gridpoints, size_t gridsize, size_t dimensions,
                                    double *alpha, base::DataMatrix &data,
                                    double treshold, size_t k,
                                    sgpp::base::OCLOperationConfiguration *parameters,
                                    size_t platformid, size_t deviceid) {
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;
  manager = std::make_shared<base::OCLManagerMultiPlatform>((*parameters)["VERBOSE"].getBool());
  DensityOCLMultiPlatform::OperationPruneGraphOCL::load_default_parameters(parameters);

  if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
    return new DensityOCLMultiPlatform::
        OperationPruneGraphOCLMultiPlatform<float>(gridpoints, gridsize, dimensions, alpha, data,
                                                   manager, parameters,
                                                   static_cast<float>(treshold), k,
                                                   platformid, deviceid);
  } else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
    return new DensityOCLMultiPlatform::
        OperationPruneGraphOCLMultiPlatform<double>(gridpoints, gridsize, dimensions, alpha, data,
                                                    manager, parameters, treshold, k, platformid,
                                                    deviceid);
  } else {
    std::stringstream errorString;
    errorString << "Error creating operation\"OperationPruneGraphOCL\": "
                << " invalid value for parameter \"INTERNAL_PRECISION\"";
    throw base::factory_exception(errorString.str().c_str());
  }
  return nullptr;
}
DensityOCLMultiPlatform::OperationPruneGraphOCL*
pruneNearestNeighborGraphConfigured(base::Grid& grid, size_t dimensions, base::DataVector &alpha,
                                    base::DataMatrix &data, double treshold, size_t k,
                                    std::string opencl_conf) {
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;

  std::cout << "Using configuration file " << opencl_conf << std::endl;
  sgpp::base::OCLOperationConfiguration *parameters =
      new sgpp::base::OCLOperationConfiguration(opencl_conf);
  manager = std::make_shared<base::OCLManagerMultiPlatform>((*parameters)["VERBOSE"].getBool());
  DensityOCLMultiPlatform::OperationPruneGraphOCL::load_default_parameters(parameters);

  size_t platformid = 0;
  if (parameters->contains("USE_PLATFORM") == true) {
    platformid = (*parameters)["USE_PLATFORM"].getInt();
  } else {
    std::stringstream errorString;
    errorString << "Error creating operation\"OperationPruneGraphOCL\": "
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
    errorString << "Error creating operation\"OperationPruneGraphOCL\": "
                << "There is no given information about which opencl device (ID) to use!"
                << " Add \"USE_DEVICE\": device_id to your configuration file "
                << "or use a different factory method." << std::endl;
    throw base::factory_exception(errorString.str().c_str());
  }

  if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
    return new DensityOCLMultiPlatform::
        OperationPruneGraphOCLMultiPlatform<float>(grid, alpha, data, dimensions, manager,
                                                   parameters,
                                                   static_cast<float>(treshold), k,
                                                   platformid, deviceid);
  } else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
    return new DensityOCLMultiPlatform::
        OperationPruneGraphOCLMultiPlatform<double>(grid, alpha, data, dimensions, manager,
                                                    parameters, treshold, k,
                                                    platformid, deviceid);
  } else {
    std::stringstream errorString;
    errorString << "Error creating operation\"OperationPruneGraphOCL\": "
                << " invalid value for parameter \"INTERNAL_PRECISION\"";
    throw base::factory_exception(errorString.str().c_str());
  }
  return nullptr;
}
}  // namespace datadriven
}  // namespace sgpp
