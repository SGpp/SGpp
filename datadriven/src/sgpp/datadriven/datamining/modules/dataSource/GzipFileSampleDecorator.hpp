// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef ZLIB
#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleDecorator.hpp>

#include <string>
#include <vector>

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
   * Decompresses a .gz file and delegates the contents down to the
   * sample provider.
   * @param fileName path to the file
   * @param hasTargets whether the file has targets (i.e. supervised learning)
   * @param readinCutoff see FileSampleProvider.hpp
   * @param readinColumns see FileSampleProvider.hpp
   * @param readinClasses see FileSampleProvider.hpp
   */
  void readFile(const std::string &fileName,
                bool hasTargets,
                size_t readinCutoff = -1,
                std::vector<size_t> readinColumns = std::vector<size_t>(),
                std::vector<double> readinClasses = std::vector<double>()) override;

  /**
   * Resets the state of the sample provider (e.g. to start a new epoch)
   */
  void reset() override;
};

} /* namespace datadriven */
} /* namespace sgpp */
#endif
