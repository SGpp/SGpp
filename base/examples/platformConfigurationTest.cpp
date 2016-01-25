/*
 * platformConfigurationTest.cpp
 *
 *  Created on: Nov 17, 2015
 *      Author: pfandedd
 */

#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

int main(int argc, char **argv) {
//    SGPP::base::OCLManagerMultiPlatform manager;

//    std::shared_ptr<SGPP::base::OCLOperationConfiguration> configuration = manager.getConfiguration();
//
//    configuration->serialize("detectPlatform.cfg");

    std::shared_ptr<SGPP::base::OCLOperationConfiguration> configuration = std::make_shared<SGPP::base::OCLOperationConfiguration>("detectPlatform.cfg");
    (*configuration).replaceIDAttr("VERBOSE", true);

    SGPP::base::OCLManagerMultiPlatform manager(configuration);

    configuration->serialize("detectPlatformOut.cfg");



}
