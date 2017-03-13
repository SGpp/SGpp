#pragma once

#include <hpx/include/components.hpp>

namespace sgpp {
namespace datadriven {
namespace MultipleEvalHPX {

struct LoadBalancerComponent
    : hpx::components::component_base<sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent> {
  LoadBalancerComponent() : result(0) {}

  LoadBalancerComponent(int32_t, size_t resultSize) : result(resultSize) {}

  std::vector<size_t> get_work_segment(size_t scheduleSize, size_t blockSize) { return {0, 0, 0}; };

  HPX_DEFINE_COMPONENT_ACTION(sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent,
                              get_work_segment, get_work_segment_action);

  void send_result_segment(std::string resultString, size_t startIndexData, size_t endIndexData){

  };

  HPX_DEFINE_COMPONENT_ACTION(sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent,
                              send_result_segment, send_result_segment_action);

  std::vector<double> get_result() { return result; }

  HPX_DEFINE_COMPONENT_ACTION(sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent, get_result,
                              get_result_action);

 private:
  std::vector<double> result;
};
}
}
}

HPX_REGISTER_ACTION_DECLARATION(
    sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent::get_work_segment_action);
