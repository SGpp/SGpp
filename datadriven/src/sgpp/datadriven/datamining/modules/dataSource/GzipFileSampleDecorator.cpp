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

#ifdef ZLIB

#include <sgpp/datadriven/datamining/modules/dataSource/GzipFileSampleDecorator.hpp>

#include <sgpp/base/exception/file_exception.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/SampleProvider.hpp>

#include <zlib.h>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

GzipFileSampleDecorator::GzipFileSampleDecorator(FileSampleProvider* const fileSampleProvider)
    : FileSampleDecorator(fileSampleProvider) {}

SampleProvider* GzipFileSampleDecorator::clone() const {
  return dynamic_cast<SampleProvider*>(new GzipFileSampleDecorator{*this});
}

void GzipFileSampleDecorator::readFile(const std::string& fileName) {
  gzFile inFileZ = gzopen(fileName.c_str(), "rb");

  if (inFileZ == nullptr) {
    throw base::file_exception("failed to open Gzip compressed file.");
  }

  size_t unzippedBytes = 0;
  std::vector<char> unzippedData(8192);

  std::string convert;
  convert.reserve(unzippedData.size());

  while (true) {
    unzippedBytes =
        gzread(inFileZ, unzippedData.data(), static_cast<unsigned int>(unzippedData.size() - 1));

    if (unzippedBytes > 0) {
      auto& last = convert.back();
      last = '\0';
      convert.append(unzippedData.data());

    } else {
      break;
    }
  }

  gzclose(inFileZ);

  fileSampleProvider->readString(convert, false);
}

void GzipFileSampleDecorator::reset() {
}

} /* namespace datadriven */
} /* namespace sgpp */
#endif
