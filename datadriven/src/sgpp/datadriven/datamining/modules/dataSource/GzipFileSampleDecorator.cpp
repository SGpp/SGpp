/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * GzipFileSampleDecorator.cpp
 *
 *  Created on: 01.04.2016
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/modules/dataSource/GzipFileSampleDecorator.hpp>

#include <sgpp/base/exception/file_exception.hpp>

#include <zlib.h>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

GzipFileSampleDecorator::GzipFileSampleDecorator(FileSampleProvider* fileSampleProvider)
    : FileSampleDecorator(fileSampleProvider) {}

GzipFileSampleDecorator::~GzipFileSampleDecorator() {}

Dataset* GzipFileSampleDecorator::getNextSamples(size_t howMany) {
  return fileSampleProvider->getNextSamples(howMany);
}

Dataset* GzipFileSampleDecorator::getAllSamples() { return fileSampleProvider->getAllSamples(); }

size_t GzipFileSampleDecorator::getDim() { return fileSampleProvider->getDim(); }

size_t GzipFileSampleDecorator::getDatasetSize() { return fileSampleProvider->getDatasetSize(); }

void GzipFileSampleDecorator::readFile(const std::string& fileName) {
  // TODO(lettrich): check if this is efficient and maybe optimize
  gzFile inFileZ = gzopen(fileName.c_str(), "rb");

  if (inFileZ == nullptr) {
    throw base::file_exception("failed to open Gzip compressed file.");
    exit(-1);
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

  fileSampleProvider->readString(convert.str());
}

void GzipFileSampleDecorator::readString(const std::string& input) {
  // TODO(lettrich): implement decompression from string, eg use
  // https://panthema.net/2007/0328-ZLibString.html as template
  fileSampleProvider->readString(input);
}

} /* namespace datadriven */
} /* namespace sgpp */
