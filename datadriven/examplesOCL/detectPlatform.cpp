// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

#include <iostream>


int main(int argc, char **argv) {
  sgpp::base::OCLManagerMultiPlatform manager(true);

  auto configuration = manager.getConfiguration();

  configuration->serialize("detectedPlatform.cfg");

  std::cout << "done" << std::endl;

  return 0;
}
