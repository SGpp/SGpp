// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/OperationMultiEvalHPX/LoadBalancerComponent.hpp>

#include <hpx/include/iostreams.hpp>
#include <hpx/include/lcos.hpp>
#include <iostream>

HPX_REGISTER_COMPONENT(
    hpx::components::component<sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent>,
    LoadBalancerComponent);

HPX_REGISTER_ACTION(
    sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent::get_work_segment_action);

HPX_REGISTER_ACTION(
    sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent::send_result_segment_action);

HPX_REGISTER_ACTION(sgpp::datadriven::MultipleEvalHPX::LoadBalancerComponent::get_result_action);
