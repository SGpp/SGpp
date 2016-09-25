#pragma once

#include <cinttypes>
#include <sstream>
#include <hpx/include/components.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/async.hpp>

#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/base/datatypes/DataVector.hpp"
#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp"
#include "sgpp/base/exception/not_implemented_exception.hpp"
#include "LoadBalancerComponent.hpp"

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
    uint32_t managerId;

    // TODO: why does this get called?
    LocalityMultiplier() :
            managerId(0) {
    }

    LocalityMultiplier(std::string serializedGrid,
            std::string serializedDataset, std::string parametersString,
            sgpp::datadriven::OperationMultipleEvalType type,
            sgpp::datadriven::OperationMultipleEvalSubType subType,
            uint32_t managerId) :
            managerId(managerId) {

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

    void mult(std::string alphaSerialized) {
        hpx::cout << "processing mult_fragment!" << std::endl << hpx::flush;
        sgpp::base::DataVector alpha = sgpp::base::DataVector::fromString(
                alphaSerialized);
        hpx::components::client<
                sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent> loadBalancer;

        loadBalancer.connect_to("manager#" + std::to_string(managerId));

        std::vector<hpx::future<void>> resultFutures;

        while (true) {
            std::vector<size_t> s =
                    hpx::async<
                            sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent::get_work_segment_action>(
                            loadBalancer.get_id(), 10000, 1).get();
            bool segmentAvailable = s[0];
            size_t startIndexData = s[1];
            size_t endIndexData = s[2];
            if (!segmentAvailable) {
                break;
            }
//            hpx::cout << "processing from " << startIndexData << " to " << endIndexData << " on locality " << hpx::get_locality_id() << std::endl << hpx::flush;

            sgpp::base::DataVector result;
            nodeMultiEval->mult(alpha, result, startIndexData, endIndexData);

            std::string resultString = result.toString();

            hpx::future<void> f =
                    hpx::async<
                            sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent::send_result_segment_action>(
                            loadBalancer.get_id(), resultString, startIndexData, endIndexData);
            resultFutures.push_back(std::move(f));
        }
        hpx::wait_all(resultFutures);
    }

    //TODO: add update grid action!

    HPX_DEFINE_COMPONENT_ACTION(sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier, mult,
            mult_action);

};

}
}
}

HPX_REGISTER_ACTION_DECLARATION(
        sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier::mult_action);

