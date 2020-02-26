// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_detectPlatformConfiguration_cpp Detect the configuration of OpenCL platforms
 *
 * This code detects the configuration of the OpenCL platforms available on the
 * machine and outputs it to a file.
 * It can be used as a quick check whether OpenCL is working.
 */

#if USE_OCL == 1

#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

int main(int argc, char** argv) {
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> configuration =
      std::make_shared<sgpp::base::OCLOperationConfiguration>(
          "detectPlatform.cfg");
  (*configuration).replaceIDAttr("VERBOSE", true);

  sgpp::base::OCLManagerMultiPlatform manager(configuration);

  configuration->serialize("detectPlatformOut.cfg");

  return 0;
}
#else
#include <iostream>
int main(int argc, char** argv) {
  std::cout << "no OpenCL support" << std::endl;
}
#endif
