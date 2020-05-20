// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef ZLIB
#include "testsCommon.hpp"

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

#include <zlib.h>
#include <boost/lexical_cast.hpp>

#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

std::string uncompressFile(std::string fileName) {
  gzFile inFileZ = gzopen(fileName.c_str(), "rb");

  if (inFileZ == nullptr) {
    std::cout << "Error: Failed to gzopen file " << fileName << std::endl;
    exit(0);
  }

  unsigned char unzipBuffer[8192];
  unsigned int unzippedBytes;
  std::vector<unsigned char> unzippedData;

  while (true) {
    unzippedBytes = gzread(inFileZ, unzipBuffer, 8192);

    if (unzippedBytes > 0) {
      for (size_t i = 0; i < unzippedBytes; i++) {
        unzippedData.push_back(unzipBuffer[i]);
      }
    } else {
      break;
    }
  }

  gzclose(inFileZ);

  std::stringstream convert;

  for (size_t i = 0; i < unzippedData.size(); i++) {
    convert << unzippedData[i];
  }

  return convert.str();
}

void doRandomRefinements(sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                         sgpp::base::Grid& grid, sgpp::base::GridGenerator& gridGen) {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1, 100);

  sgpp::base::DataVector alphaRefine(grid.getSize());

  for (size_t i = 0; i < alphaRefine.getSize(); i++) {
    alphaRefine[i] = dist(mt);
  }

  for (size_t i = 0; i < adaptivityConfig.numRefinements_; i++) {
    sgpp::base::SurplusRefinementFunctor myRefineFunc(
        alphaRefine, adaptivityConfig.numRefinementPoints_, adaptivityConfig.refinementThreshold_);
    gridGen.refine(myRefineFunc);
    size_t oldSize = alphaRefine.getSize();
    alphaRefine.resize(grid.getSize());

    for (size_t j = oldSize; j < alphaRefine.getSize(); j++) {
      alphaRefine[j] = dist(mt);
    }
  }
}

void doDirectedRefinements(sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                           sgpp::base::Grid& grid, sgpp::base::GridGenerator& gridGen) {
  double dummySurplusValue = 1.0;

  sgpp::base::DataVector alphaRefine(grid.getSize());

  for (size_t i = 0; i < alphaRefine.getSize(); i++) {
    alphaRefine[i] = dummySurplusValue;
    dummySurplusValue += 1.0;
  }

  for (size_t i = 0; i < adaptivityConfig.numRefinements_; i++) {
    sgpp::base::SurplusRefinementFunctor myRefineFunc(
        alphaRefine, adaptivityConfig.numRefinementPoints_, adaptivityConfig.refinementThreshold_);
    gridGen.refine(myRefineFunc);
    size_t oldSize = alphaRefine.getSize();
    alphaRefine.resize(grid.getSize());

    for (size_t j = oldSize; j < alphaRefine.getSize(); j++) {
      alphaRefine[i] = dummySurplusValue;
    }

    // increment only once for added grid points
    dummySurplusValue += 1.0;
  }
}

#endif
