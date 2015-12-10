/*
 * kernelTuner.cpp
 *
 *  Created on: Oct 27, 2015
 *      Author: pfandedd
 */
#include <iostream>

#if USE_OCL == 1

#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>


int main(int argc, char** argv) {

  SGPP::base::OCLManagerMultiPlatform manager;

  auto configuration = manager.getConfiguration();

  configuration->serialize("detectedPlatform.cfg");

  std::cout << "done" << std::endl;

  return 0;
}
#else
int main(int argc, char** argv) {
  std::cout << "no OpenCL support" << std::endl;
}
#endif

