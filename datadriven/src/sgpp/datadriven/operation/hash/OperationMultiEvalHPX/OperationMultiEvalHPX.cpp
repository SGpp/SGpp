// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <hpx/include/async.hpp>
#include <hpx/include/components.hpp>
#include <hpx/include/lcos.hpp>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include <sgpp/datadriven/operation/hash/OperationMultiEvalHPX/LocalityMultiplier.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultiEvalHPX/OperationMultiEvalHPX.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/base/tools/QueueLoadBalancerMutex.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultiEvalStreaming/OperationMultiEvalStreaming.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/OperatorFactory.hpp>

namespace sgpp {
namespace datadriven {

OperationMultiEvalHPX::OperationMultiEvalHPX(
    base::Grid& grid, base::DataMatrix& dataset,
    sgpp::datadriven::OperationMultipleEvalConfiguration& configuration, bool verbose)
    : OperationMultipleEval(grid, dataset),
      configuration(configuration),
      dim(grid.getDimension()),
      verbose(verbose),
      duration(-1.0) {
  // create the kernel specific data structures for the current grid
  this->prepare();
}

OperationMultiEvalHPX::~OperationMultiEvalHPX() {}

void OperationMultiEvalHPX::mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  this->myTimer.start();

  result.resize(dataset.getNrows());

  std::vector<hpx::id_type> all_ids = hpx::find_remote_localities();

  hpx::default_distribution_policy policy = hpx::default_layout(all_ids);

  std::string serializedGrid;
  grid.serialize(serializedGrid);

  std::string serializedDataset;
  dataset.toString(serializedDataset);

  std::ostringstream transferStream;
  configuration.getParameters()->serialize(transferStream, 0);
  std::string transferableParameterString = transferStream.str();

  hpx::components::client<sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent> loadBalancer =
      hpx::new_<hpx::components::client<sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent>>(
          hpx::find_here(), 0, dataset.getNrows());
  loadBalancer.register_as("manager#" + std::to_string(hpx::get_locality_id()), false);

  std::vector<hpx::components::client<sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier>>
      multipliers =
          hpx::new_<
              hpx::components::client<sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier>[]>(
              policy, all_ids.size(), serializedGrid, serializedDataset,
              transferableParameterString, configuration.getType(), configuration.getSubType(),
              hpx::get_locality_id())
              .get();

  std::string alphaSerialized = alpha.toString();

  std::vector<hpx::future<void>> finished;
  for (hpx::components::client<sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier>& multiplier :
       multipliers) {
    hpx::future<void> f =
        hpx::async<sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier::mult_action>(
            multiplier.get_id(), alphaSerialized);
    finished.push_back(std::move(f));
  }

  hpx::wait_all(finished);
  std::vector<double> resultVector =
      hpx::async<sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent::get_result_action>(
          loadBalancer.get_id())
          .get();
  for (size_t i = 0; i < resultVector.size(); i++) {
    result[i] = static_cast<double>(resultVector.at(i));
  }

  this->duration = this->myTimer.stop();
}

void OperationMultiEvalHPX::multTranspose(sgpp::base::DataVector& source,
                                          sgpp::base::DataVector& result) {
  this->myTimer.start();
  this->duration = this->myTimer.stop();
  throw base::not_implemented_exception();
}

double OperationMultiEvalHPX::getDuration() { return this->duration; }

void OperationMultiEvalHPX::prepare() {}
}  // namespace datadriven
}  // namespace sgpp
