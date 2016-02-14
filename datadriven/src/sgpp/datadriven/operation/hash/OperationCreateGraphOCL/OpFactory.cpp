#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/Configuration.hpp>
#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>
#include "OpFactory.hpp"
namespace SGPP {
namespace datadriven {

SGPP::datadriven::StreamingOCLMultiPlatform::OperationCreateGraphOCL* createNearestNeighbarGraphConfigured(size_t k, size_t dimensions, std::string opencl_conf) {
	std::shared_ptr<base::OCLManagerMultiPlatform> manager;

	std::cout<<"Using configuration file "<<opencl_conf<<std::endl;
	SGPP::base::OCLOperationConfiguration *parameters=new SGPP::base::OCLOperationConfiguration(opencl_conf);
	manager = std::make_shared<base::OCLManagerMultiPlatform>(true);
	parameters->serialize("MyOCLConfDebug.cfg");
	if(parameters->contains("INTERNAL_PRECISION")==false) {
		std::cout<<"Warning! No internal precision setting detected. Using single precision from now on!"<<std::endl;
		parameters->addIDAttr("INTERNAL_PRECISION", "float");
	}
	if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
	}else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
	} else {
		throw base::factory_exception(
			"Error creating operation\"CreateGraphOCL\": invalid value for parameter \"INTERNAL_PRECISION\"");
	}
	parameters->serialize("MyOCLConf.cfg");

	std::string &firstPlatformName = (*parameters)["PLATFORMS"].keys()[0];
	std::string &firstDeviceName = (*parameters)["PLATFORMS"][firstPlatformName]["DEVICES"].keys()[0];
	json::Node &deviceNode = (*parameters)["PLATFORMS"][firstPlatformName]["DEVICES"][firstDeviceName];
	json::Node &firstDeviceConfig = deviceNode["KERNELS"]["multdensity"];

	if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
return new SGPP::datadriven::StreamingOCLMultiPlatform::OperationCreateGraphOCLMultiPlatform<float>(dimensions, manager, deviceNode, k);
	} else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
return new SGPP::datadriven::StreamingOCLMultiPlatform::OperationCreateGraphOCLMultiPlatform<double>(dimensions, manager, deviceNode, k);
	} else {
		throw base::factory_exception(
			"Error creating operation\"CreateGraphOCL\": invalid value for parameter \"INTERNAL_PRECISION\"");
	}
	return NULL;
}
}
}
