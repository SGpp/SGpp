// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <algorithm>
#include <iostream>
#include <vector>
#include <sstream>

#include "OperationMultiEvalHPX.hpp"
#include "LocalityMultiplier.hpp"
#include "sgpp/base/exception/not_implemented_exception.hpp"
#include "sgpp/datadriven/operation/hash/OperationMultiEvalStreaming/OperationMultiEvalStreaming.hpp"
#include "sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/OperatorFactory.hpp"
#include "sgpp/base/opencl/QueueLoadBalancer.hpp"

#include <hpx/include/components.hpp>
#include <hpx/include/async.hpp>

namespace sgpp {
namespace datadriven {

OperationMultiEvalHPX::OperationMultiEvalHPX(base::Grid& grid,
        base::DataMatrix& dataset,
        sgpp::datadriven::OperationMultipleEvalConfiguration& configuration,
        bool verbose) :
        OperationMultipleEval(grid, dataset), configuration(configuration), dim(
                grid.getDimension()), verbose(verbose), duration(-1.0) {
    // create the kernel specific data structures for the current grid
    this->prepare();
}

OperationMultiEvalHPX::~OperationMultiEvalHPX() {
}

void OperationMultiEvalHPX::mult(sgpp::base::DataVector& alpha,
        sgpp::base::DataVector& result) {
    this->myTimer.start();

    result.resize(dataset.getNrows());

//    // create appropriate node level multi eval implementation
//    std::unique_ptr<sgpp::base::OperationMultipleEval> nodeMultiEval;
//    if (configuration.getType() == OperationMultipleEvalType::STREAMING
//            && configuration.getSubType()
//                    == OperationMultipleEvalSubType::DEFAULT) {
//        nodeMultiEval =
//                std::make_unique<datadriven::OperationMultiEvalStreaming>(grid,
//                        dataset);
//    } else if (configuration.getType() == OperationMultipleEvalType::STREAMING
//            && configuration.getSubType()
//                    == OperationMultipleEvalSubType::OCLMP) {
//        nodeMultiEval = std::unique_ptr<OperationMultipleEval>(
//                createStreamingOCLMultiPlatformConfigured(grid, dataset,
//                        configuration));
//    } else {
//        throw base::not_implemented_exception();
//    }
//    nodeMultiEval->mult(alpha, result);

    size_t num_localities = hpx::get_num_localities().get();
    std::vector<hpx::id_type> all_ids = hpx::find_remote_localities();

    hpx::default_distribution_policy policy = hpx::default_layout(all_ids);

    std::string serializedGrid;
    grid.serialize(serializedGrid);

    std::string serializedDataset;
    dataset.toString(serializedDataset);

    std::ostringstream transferStream;
    configuration.getParameters()->serialize(transferStream, 0);
    std::string transferableParameterString = transferStream.str();

//    std::vector<
//            hpx::components::client<
//                    sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier<float>>>multipliers =
//    hpx::new_<hpx::components::client<sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier<float>>[]>(
//            policy, num_localities, serializedGrid, serializedDataset, transferableParameterString, configuration.getType(), configuration.getSubType()).get();
//
//    for (hpx::components::client<
//            sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier<float>> &multiplier : multipliers) {
//        uint32_t comp_locality = hpx::naming::get_locality_id_from_id(
//                multiplier.get_id());
//        multiplier.register_as("/multiplier#" + std::to_string(comp_locality),
//                false);
//    }

    hpx::components::client<
            sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier> multiplier =
            hpx::new_<
                    hpx::components::client<
                            sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier>>(
                    policy, serializedGrid, serializedDataset,
                    transferableParameterString, configuration.getType(),
                    configuration.getSubType());
    multiplier.wait();
    std::string alphaSerialized = alpha.toString();

    base::QueueLoadBalancer queueLoadBalancer;
    queueLoadBalancer.initialize(0, dataset.getNrows());

    size_t startIndex = 0;
    size_t stopIndex = 0;
    const size_t blockSize = 1;
    const size_t scheduleSize = 10000;
    while (true) {
        bool segmentAvailable = queueLoadBalancer.getNextSegment(scheduleSize,
                blockSize, startIndex, stopIndex);
        if (!segmentAvailable) {
            break;
        }
        hpx::future<std::string> f =
                hpx::async<
                        sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier::mult_fragment_action>(
                        multiplier.get_id(), alphaSerialized, startIndex,
                        stopIndex);
        f.wait();
        hpx::future<void> g = f.then(
                hpx::util::unwrapped(
                        [=,&result](std::string resultSerialized)
                        {
                            sgpp::base::DataVector resultSegment = sgpp::base::DataVector::fromString(resultSerialized);
                            for (size_t i = 0; i < resultSegment.getSize(); i++) {
                                result[startIndex + i] = resultSegment[i];
                            }
                        }));
        g.wait();
    }

    this->duration = this->myTimer.stop();
}

void OperationMultiEvalHPX::multTranspose(sgpp::base::DataVector& source,
        sgpp::base::DataVector& result) {
    this->myTimer.start();
    this->duration = this->myTimer.stop();
    throw base::not_implemented_exception();
}

double OperationMultiEvalHPX::getDuration() {
    return this->duration;
}

void OperationMultiEvalHPX::prepare() {
}
}  // namespace datadriven
}  // namespace sgpp
