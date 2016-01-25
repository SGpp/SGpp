/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalAdaptiveOCL/OperationMultiEvalAdaptiveOCL.hpp>

namespace SGPP {
  namespace datadriven {

    base::OperationMultipleEval* createAdaptiveOCLConfigured(base::Grid& grid, base::DataMatrix& dataset,
        SGPP::datadriven::OperationMultipleEvalConfiguration& configuration) {

      std::shared_ptr<base::OCLOperationConfiguration> parameters;

      if (configuration.getParameters().operator bool()) {
        base::OCLOperationConfiguration* cloned =
          dynamic_cast<base::OCLOperationConfiguration*>(configuration.getParameters()->clone());
        parameters = std::shared_ptr<base::OCLOperationConfiguration>(cloned);
      } else {
        parameters = std::make_shared<base::OCLOperationConfiguration>("AdaptiveOCL.cfg");

        if ((*parameters).contains("KERNEL_USE_LOCAL_MEMORY") == false) {
          (*parameters).addIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
        }

        if ((*parameters).contains("LOCAL_SIZE") == false) {
          (*parameters).addIDAttr("LOCAL_SIZE", 128ul);
        }

        if ((*parameters).contains("KERNEL_DATA_BLOCKING_SIZE") == false) {
          (*parameters).addIDAttr("KERNEL_DATA_BLOCKING_SIZE", 1ul);
        }

        if ((*parameters).contains("LINEAR_LOAD_BALANCING_VERBOSE") == false) {
          (*parameters).addIDAttr("LINEAR_LOAD_BALANCING_VERBOSE", false);
        }

        if ((*parameters).contains("KERNEL_TRANS_DATA_BLOCK_SIZE") == false) {
          (*parameters).addIDAttr("KERNEL_TRANS_DATA_BLOCK_SIZE", 1ul);
        }

        if ((*parameters).contains("ADAPTIVE_STREAMING_HARD_LIMIT") == false) {
          (*parameters).addIDAttr("ADAPTIVE_STREAMING_HARD_LIMIT", 10ul); //absolute value
        }

        if ((*parameters).contains("ADAPTIVE_STREAMING_DENSITY") == false) {
          (*parameters).addIDAttr("ADAPTIVE_STREAMING_DENSITY", 5ul); //In percent
        }
      }

      if ((*parameters)["VERBOSE"].getBool()) {
        std::cout << "are optimizations on: " << (*parameters)["ENABLE_OPTIMIZATIONS"].getBool() << std::endl;
        std::cout << "is local memory on: " << (*parameters)["KERNEL_USE_LOCAL_MEMORY"].getBool() << std::endl;
        std::cout << "local size: " << (*parameters)["LOCAL_SIZE"].getUInt() << std::endl;
        std::cout << "internal precision: " << (*parameters)["INTERNAL_PRECISION"].get() << std::endl;
        std::cout << "platform is: " << (*parameters)["PLATFORM"].get() << std::endl;
        std::cout << "device type is: " << (*parameters)["DEVICE_TYPE"].get() << std::endl;
        std::cout << "hard limit: " << (*parameters)["ADAPTIVE_STREAMING_HARD_LIMIT"].get() << std::endl;
        std::cout << "soft limit: " << (*parameters)["ADAPTIVE_STREAMING_DENSITY"].get() << std::endl;
      }

      if ((*parameters)["INTERNAL_PRECISION"].get() == "float") {
        return new datadriven::OperationMultiEvalAdaptiveOCL<float>(grid, dataset, parameters);
      } else if ((*parameters)["INTERNAL_PRECISION"].get() == "double") {
        return new datadriven::OperationMultiEvalAdaptiveOCL<double>(grid, dataset, parameters);
      } else {
        throw base::factory_exception(
          "Error creating operation\"OperationMultiEvalAdaptiveOCL\": invalid value for parameter \"INTERNAL_PRECISION\"");
      }

      //TODO: make parameters changable through api
      /*std::shared_ptr<base::OCLOperationConfiguration> parameters;
       (*parameters).addTextAttr("KERNEL_USE_LOCAL_MEMORY", "true");
       (*parameters).addTextAttr("KERNEL_MAX_DIM_UNROLL", "10");
       (*parameters).addTextAttr("LINEAR_LOAD_BALANCING_VERBOSE", "false");

       parameters->readFromFile("AdaptiveOCL.cfg");

       //  std::cout << "are optimizations on: " << (*parameters)["STREAMING_OCL_ENABLE_OPTIMIZATIONS"].getBool() << std::endl;
       //  std::cout << "is local memory on: " << (*parameters)["STREAMING_OCL_USE_LOCAL_MEMORY"].getBool() << std::endl;
       //  std::cout << "local size: " << (*parameters)["STREAMING_OCL_LOCAL_SIZE"].getUInt() << std::endl;
       //  std::cout << "max dim unroll: " << (*parameters)["STREAMING_OCL_MAX_DIM_UNROLL"].getUInt() << std::endl;
       //  std::cout << "internal precision: " << (*parameters)["STREAMING_OCL_INTERNAL_PRECISION"].get() << std::endl;
       //  std::cout << "platform is: " << (*parameters)["STREAMING_OCL_PLATFORM"].get() << std::endl;
       //  std::cout << "device type is: " << (*parameters)["STREAMING_OCL_DEVICE_TYPE"].get() << std::endl;

       if ((*parameters)["INTERNAL_PRECISION"].get() == "float") {
       return new datadriven::OperationMultiEvalAdaptiveOCL<float>(grid, dataset, parameters);
       } else if ((*parameters)["INTERNAL_PRECISION"].get() == "double") {
       return new datadriven::OperationMultiEvalAdaptiveOCL<double>(grid, dataset, parameters);
       } else {
       throw base::factory_exception(
       "Error creating operation\"OperationMultiEvalAdaptiveOCL\": invalid value for parameter \"INTERNAL_PRECISION\"");
       }*/
    }

  }
}
