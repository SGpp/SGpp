#pragma once

#include <cinttypes>
#include <sstream>
#include <hpx/include/components.hpp>
#include <hpx/include/iostreams.hpp>

#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/base/datatypes/DataVector.hpp"
#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp"
#include "sgpp/base/exception/not_implemented_exception.hpp"

#include "sgpp/datadriven/operation/hash/OperationMultiEvalStreaming/OperationMultiEvalStreaming.hpp"
#include "sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/OperatorFactory.hpp"

namespace sgpp {
namespace datadriven {
namespace MultipleEvalHPX {

struct LocalityMultiplier: hpx::components::component_base<
        sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier> {
    std::unique_ptr<sgpp::base::Grid> grid;
    sgpp::base::DataMatrix dataset;
    sgpp::datadriven::OperationMultipleEvalConfiguration configuration;
    std::unique_ptr<sgpp::base::OperationMultipleEval> nodeMultiEval;

    // TODO: why does this get called?
    LocalityMultiplier() {
    }

    LocalityMultiplier(std::string serializedGrid,
            std::string serializedDataset, std::string parametersString,
            sgpp::datadriven::OperationMultipleEvalType type,
            sgpp::datadriven::OperationMultipleEvalSubType subType) {

        std::unique_ptr<sgpp::base::OCLOperationConfiguration> parameters =
                std::make_unique<sgpp::base::OCLOperationConfiguration>();
        hpx::cout << parametersString << std::endl << hpx::flush;
        parameters->deserializeFromString(parametersString);

        grid = std::unique_ptr<sgpp::base::Grid>(
                sgpp::base::Grid::unserialize(serializedGrid));
        dataset = sgpp::base::DataMatrix::fromString(serializedDataset);

        //TODO: have to add toVector -> fromVector to LinearGrid

        // create node-level operation configuration
        sgpp::datadriven::OperationMultipleEvalConfiguration configuration(type,
                subType, sgpp::datadriven::OperationMultipleEvalMPIType::NONE,
                *parameters);
        nodeMultiEval = std::unique_ptr<sgpp::base::OperationMultipleEval>(
                sgpp::op_factory::createOperationMultipleEval(*grid, dataset,
                        configuration));
        hpx::cout << "component set up!" << std::endl << hpx::flush;
    }

    std::string mult_fragment(std::string alphaSerialized, size_t startIndexData, size_t endIndexData) {
        hpx::cout << "processing mult_fragment!" << std::endl << hpx::flush;
        sgpp::base::DataVector alpha = sgpp::base::DataVector::fromString(
                alphaSerialized);

        sgpp::base::DataVector result;
//        result.resize(endIndexData - startIndexData);
//        result.setAll(0.0);
        nodeMultiEval->mult(alpha, result, startIndexData, endIndexData);

        std::string resultString = result.toString();
        return resultString;
    }

    //TODO: add update grid action!

    HPX_DEFINE_COMPONENT_ACTION(sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier, mult_fragment,
            mult_fragment_action);

};

}
}
}

HPX_REGISTER_ACTION_DECLARATION(
        sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier::mult_fragment_action);

