#include "LocalityMultiplier.hpp"

#include <iostream>
#include <hpx/include/lcos.hpp>
#include <hpx/include/iostreams.hpp>

HPX_REGISTER_COMPONENT(
        hpx::components::component<
                sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier>,
        LocalityMultiplier);

HPX_REGISTER_ACTION(
        sgpp::datadriven::MultipleEvalHPX::LocalityMultiplier::mult_fragment_action);
