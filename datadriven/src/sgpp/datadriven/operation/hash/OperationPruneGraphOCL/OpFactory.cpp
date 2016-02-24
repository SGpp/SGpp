// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/globaldef.hpp>
#include <string>
#include "OpFactory.hpp"
#include "KernelPruneGraph.hpp"
namespace SGPP {
namespace datadriven {

DensityOCLMultiPlatform::OperationPruneGraphOCL*
pruneNearestNeighborGraphConfigured(base::Grid& grid, size_t dimensions, base::DataVector &alpha,
                                    base::DataMatrix &data, double treshold, size_t k,
                                    std::string opencl_conf) {
    std::shared_ptr<base::OCLManagerMultiPlatform> manager;

    std::cout << "Using configuration file " << opencl_conf << std::endl;
    SGPP::base::OCLOperationConfiguration *parameters =
        new SGPP::base::OCLOperationConfiguration(opencl_conf);
    manager = std::make_shared<base::OCLManagerMultiPlatform>(true);
    parameters->serialize("MyOCLConfDebug.cfg");
    if (parameters->contains("INTERNAL_PRECISION") == false) {
        std::cout << "Warning! No internal precision setting detected."
                  << " Using double precision from now on!" << std::endl;
        parameters->addIDAttr("INTERNAL_PRECISION", "float");
    }
    if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
        DensityOCLMultiPlatform::KernelPruneGraph<float>::augmentDefaultParameters(*parameters);
    } else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
        DensityOCLMultiPlatform::KernelPruneGraph<double>::augmentDefaultParameters(*parameters);
    } else {
        std::stringstream errorString;
        errorString << "Error creating operation\"OperationPruneGraphOCL\": "
                    << " invalid value for parameter \"INTERNAL_PRECISION\"";
        throw base::factory_exception(errorString.str().c_str());
    }
    parameters->serialize("MyOCLConf.cfg");

    std::string &firstPlatformName =
        (*parameters)["PLATFORMS"].keys()[0];
    std::string &firstDeviceName =
        (*parameters)["PLATFORMS"][firstPlatformName]["DEVICES"].keys()[0];
    json::Node &deviceNode =
        (*parameters)["PLATFORMS"][firstPlatformName]["DEVICES"][firstDeviceName];
    json::Node &firstDeviceConfig = deviceNode["KERNELS"]["removeEdges"];

    if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
        return new DensityOCLMultiPlatform::
            OperationPruneGraphOCLMultiPlatform<float>(grid, alpha, data, dimensions, manager,
                                                       firstDeviceConfig,
                                                       static_cast<float>(treshold), k);
    } else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
        return new DensityOCLMultiPlatform::
            OperationPruneGraphOCLMultiPlatform<double>(grid, alpha, data, dimensions, manager,
                                                        firstDeviceConfig, treshold, k);
    } else {
        std::stringstream errorString;
        errorString << "Error creating operation\"OperationPruneGraphOCL\": "
                    << " invalid value for parameter \"INTERNAL_PRECISION\"";
        throw base::factory_exception(errorString.str().c_str());
    }
    return NULL;
}
}  // namespace datadriven
}  // namespace SGPP
