#pragma once

#include <sgpp/globaldef.hpp>

#include "OperationCreateGraphOCLMultiPlatform.hpp"
namespace SGPP {
namespace datadriven {

SGPP::datadriven::StreamingOCLMultiPlatform::OperationCreateGraphOCL* createNearestNeighbarGraphConfigured(SGPP::base::DataVector &data,std::string opencl_conf);
}
}
