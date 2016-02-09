#pragma once

#include <sgpp/globaldef.hpp>

#include "OperationDensityOCLMultiPlatform.hpp"
namespace SGPP {
namespace datadriven {

SGPP::datadriven::StreamingOCLMultiPlatform::OperationDensityOCL* createDensityOCLMultiPlatformConfigured(base::Grid& grid,std::string opencl_conf);
}
}
