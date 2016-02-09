/*
 * test_datadrivenCommon.cpp
 *
 *  Created on: Jul 30, 2015
 *      Author: pfandedd
 */

#include "testsCommon.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <random>
#include <fstream>

#include <boost/lexical_cast.hpp>

#include <zlib.h>

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

using namespace SGPP::base;

std::string uncompressFile(std::string fileName) {

  gzFile inFileZ = gzopen(fileName.c_str(), "rb");

  if (inFileZ == NULL) {
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

void doRandomRefinements(SGPP::base::AdpativityConfiguration& adaptConfig,
                         SGPP::base::Grid& grid, SGPP::base::GridGenerator& gridGen) {

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1, 100);

  SGPP::base::DataVector alphaRefine(grid.getSize());

  for (size_t i = 0; i < alphaRefine.getSize(); i++) {
    alphaRefine[i] = dist(mt);
  }

  for (size_t i = 0; i < adaptConfig.numRefinements_; i++) {
    SGPP::base::SurplusRefinementFunctor* myRefineFunc = new
    SGPP::base::SurplusRefinementFunctor(&alphaRefine,
                                         adaptConfig.noPoints_, adaptConfig.threshold_);
    gridGen.refine(myRefineFunc);
    size_t oldSize = alphaRefine.getSize();
    alphaRefine.resize(grid.getSize());

    for (size_t j = oldSize; j < alphaRefine.getSize(); j++) {
      alphaRefine[j] = dist(mt);
    }

    delete myRefineFunc;
  }
}

void doDirectedRefinements(SGPP::base::AdpativityConfiguration& adaptConfig,
                           SGPP::base::Grid& grid, SGPP::base::GridGenerator& gridGen) {

  double dummySurplusValue = 1.0;

  SGPP::base::DataVector alphaRefine(grid.getSize());

  for (size_t i = 0; i < alphaRefine.getSize(); i++) {
    alphaRefine[i] = dummySurplusValue;
    dummySurplusValue += 1.0;
  }

  for (size_t i = 0; i < adaptConfig.numRefinements_; i++) {
    SGPP::base::SurplusRefinementFunctor* myRefineFunc = new
    SGPP::base::SurplusRefinementFunctor(&alphaRefine,
                                         adaptConfig.noPoints_, adaptConfig.threshold_);
    gridGen.refine(myRefineFunc);
    size_t oldSize = alphaRefine.getSize();
    alphaRefine.resize(grid.getSize());

    for (size_t j = oldSize; j < alphaRefine.getSize(); j++) {
      alphaRefine[i] = dummySurplusValue;

    }

    //increment only once for added grid points
    dummySurplusValue += 1.0;

    delete myRefineFunc;
  }
}
