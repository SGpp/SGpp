#include "LoadBalancerComponent.hpp"

#include <iostream>
#include <hpx/include/lcos.hpp>
#include <hpx/include/iostreams.hpp>

HPX_REGISTER_COMPONENT(
        hpx::components::component<
                sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent>,
        LoadBalancerComponent);

HPX_REGISTER_ACTION(
        sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent::get_work_segment_action);

HPX_REGISTER_ACTION(
        sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent::send_result_segment_action);

HPX_REGISTER_ACTION(
        sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent::get_result_action);
