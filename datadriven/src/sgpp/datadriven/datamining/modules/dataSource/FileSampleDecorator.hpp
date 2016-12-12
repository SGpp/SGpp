/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FileSampleDecorator.hpp
 *
 *  Created on: 01.04.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleProvider.hpp>

namespace sgpp {
namespace datadriven {

/**
 * #FileSampleDecorator provides an interface to provide generic manipulations for various kinds of
 * #sgpp::datadriven::FileSampleProvider using the decorator pattern.
 *
 * Using inheritance and deligation we can manipulate input or output of
 * #sgpp::datadriven::FileSampleProvider member functions  in a generic fashion e.g. apply
 * decompression of compressed files before passing them to the actual implementations without
 * limiting ourselves to one specific #sgpp::datadriven::FileSampleProvider. For better
 * understanding see derived classes.
 */
class FileSampleDecorator : public FileSampleProvider {
 public:
  /**
   * Constructor.
   *
   * @param fileSampleProvider: Pointer to a sgpp::datadriven::FileSampleProvider that will be
   * wrapped as a deligate. The decorator will take ownership
   * of this object and take care of its destruction.
   */
  explicit FileSampleDecorator(FileSampleProvider* const fileSampleProvider);

  FileSampleDecorator(const FileSampleDecorator& rhs);

  FileSampleDecorator(FileSampleDecorator&& rhs) = default;

  FileSampleDecorator& operator=(const FileSampleDecorator& rhs);

  FileSampleDecorator& operator=(FileSampleDecorator&& rhs) = default;

  ~FileSampleDecorator() = default;

  Dataset* getNextSamples(size_t howMany) override;

  Dataset* getAllSamples() override;

  size_t getDim() const override;

  size_t getNumSamples() const override;

  void readFile(const std::string& fileName) override;

  void readString(const std::string& input) override;

 protected:
  /**
   * Delegate sgpp::datadriven::FileSampleProvider object. Calls to the object will be wrapped by
   * the decorator performing preprocessing or postprocessing to parts of the object, passing
   * through calls which are not of interrest.
   */
  std::unique_ptr<FileSampleProvider> fileSampleProvider;
};

} /* namespace datadriven */
} /* namespace sgpp */
