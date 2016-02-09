/*
 * platformConfigurationTest.cpp
 *
 *  Created on: Nov 17, 2015
 *      Author: pfandedd, baurms
 */
#include <iostream>

#if USE_OCL == 1

#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>


int main(int argc, char** argv) {

  std::shared_ptr<SGPP::base::OCLOperationConfiguration> configuration =
    std::make_shared<SGPP::base::OCLOperationConfiguration>("detectPlatform.cfg");
  (*configuration).replaceIDAttr("VERBOSE", true);

  SGPP::base::OCLManagerMultiPlatform manager(configuration);

  configuration->serialize("detectPlatformOut.cfg");

  return 0;
}
#else
int main(int argc, char** argv) {
  std::cout << "no OpenCL support" << std::endl;
}
#endif
