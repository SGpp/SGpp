#pragma once

#include <sgpp/globaldef.hpp>

#include "OperationDensityOCLMultiPlatform.hpp"
namespace SGPP {
namespace datadriven {

SGPP::datadriven::StreamingOCLMultiPlatform::OperationDensityOCL* createDensityOCLMultiPlatformConfigured(base::Grid& grid, size_t dimension, double lambda, std::string opencl_conf);
}
}
