// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <hpx/include/components.hpp>

#include <vector>

#include <sgpp/base/tools/QueueLoadBalancerMutex.hpp>

namespace sgpp {
namespace datadriven {
namespace MultipleEvalHPX {

struct LoadBalancerComponent
    : hpx::components::component_base<sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent> {
  // Dummy method that hpx requires (likely for the client)
  LoadBalancerComponent() : result(0) {}

  LoadBalancerComponent(int32_t, size_t resultSize) : result(resultSize) {
    balancer.initialize(0, resultSize);
  }

  std::vector<size_t> get_work_segment(size_t scheduleSize, size_t blockSize) {
    size_t segmentStart;
    size_t segmentEnd;
    bool hasNextSegment =
        balancer.getNextSegment(scheduleSize, blockSize, segmentStart, segmentEnd);
    if (hasNextSegment) {
      return {1, segmentStart, segmentEnd};
    } else {
      return {0, 0, 0};
    }
  }

  HPX_DEFINE_COMPONENT_ACTION(sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent,
                              get_work_segment, get_work_segment_action);

  void send_result_segment(std::vector<double> partialResult, size_t startIndexData,
                           size_t endIndexData) {
    for (size_t i = 0; i < endIndexData - startIndexData; i++) {
      result[i + startIndexData] = partialResult[i];
    }
  }

  HPX_DEFINE_COMPONENT_ACTION(sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent,
                              send_result_segment, send_result_segment_action);

  std::vector<double> get_result() { return result; }

  HPX_DEFINE_COMPONENT_ACTION(sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent, get_result,
                              get_result_action);

 private:
  std::vector<double> result;
  sgpp::base::QueueLoadBalancerMutex balancer;
};
}  // namespace MultipleEvalHPX
}  // namespace datadriven
}  // namespace sgpp

HPX_REGISTER_ACTION_DECLARATION(
    sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent::get_work_segment_action);
