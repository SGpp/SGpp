/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * GzipFileSampleDecorator.hpp
 *
 *  Created on: 01.04.2016
 *      Author: Michael Lettrich
 */
#ifdef ZLIB
#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleDecorator.hpp>

#include <string>

namespace sgpp {
namespace datadriven {
/**
 * Adds the ability to read gzip compressed files to file sample providers.
 *
 * This class wraps any valid #sgpp::datadriven::FileSampleProvider object and adds a decompression
 * step to the #readFile member function before trying to parse the contents of the file.
 */
class GzipFileSampleDecorator : public FileSampleDecorator {
 public:
  /**
   * Constructor decorating a FileSampleProvider object.
   *
   * @param fileSampleProvider: pointer to the object to be used as a delegate.
   */
  explicit GzipFileSampleDecorator(FileSampleProvider* fileSampleProvider);

  SampleProvider* clone() const override;

  /**
   * Decompress a gzip compressed file at the given path and read its contents using the delegate's
   * member functions. Throws if the file cannot be opened or parsed.
   *
   * @param filePath valid path to a gzip compressed file.
   */
  void readFile(const std::string& filePath);

  /**
   * Resets the state of the sample provider (e.g. to start a new epoch)
   */
  void reset() override;
};

} /* namespace datadriven */
} /* namespace sgpp */
#endif
