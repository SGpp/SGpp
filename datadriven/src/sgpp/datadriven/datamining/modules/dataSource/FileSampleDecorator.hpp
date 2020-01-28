// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleProvider.hpp>

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * #FileSampleDecorator provides an interface to provide generic manipulations for various kinds of
 * #sgpp::datadriven::FileSampleProvider using the decorator pattern.
 *
 * Using inheritance and deligation we can manipulate input or output of
 * #sgpp::datadriven::FileSampleProvider member functions  in a generic fashion e.g. apply
 * decompression of compressed files before passing them to the actual
 * #sgpp::datadriven::FileSampleProvider without limiting ourselves to a specific implementation.
 * Instead we use the decorated object's member functions as a black box.
 */
class FileSampleDecorator : public FileSampleProvider {
 public:
  /**
   * Constructor.
   *
   * @param fileSampleProvider: Pointer to a #sgpp::datadriven::FileSampleProvider that will be
   * wrapped as a deligate. The decorator will take ownership
   * of this object and take care of its destruction.
   */
  explicit FileSampleDecorator(FileSampleProvider *const fileSampleProvider);

  FileSampleDecorator(const FileSampleDecorator &rhs);

  FileSampleDecorator(FileSampleDecorator &&rhs) = default;

  FileSampleDecorator &operator=(const FileSampleDecorator &rhs);

  FileSampleDecorator &operator=(FileSampleDecorator &&rhs) = default;

  ~FileSampleDecorator() override = default;

  Dataset *getNextSamples(size_t howMany) override;

  Dataset *getAllSamples() override;

  size_t getDim() const override;

  size_t getNumSamples() const override;

  /**
   * Reads a file's content from a file
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
   * Reads a file's content
   * @param input the file's content
   * @param hasTargets whether the file has targets (i.e. supervised learning)
   * @param readinCutoff see FileSampleProvider.hpp
   * @param readinColumns see FileSampleProvider.hpp
   * @param readinClasses see FileSampleProvider.hpp
   */
  void readString(const std::string &input,
                  bool hasTargets,
                  size_t readinCutoff = -1,
                  std::vector<size_t> readinColumns = std::vector<size_t>(),
                  std::vector<double> readinClasses = std::vector<double>()) override;

 protected:
  /**
   * Delegate #sgpp::datadriven::FileSampleProvider object. Calls to the object will be wrapped by
   * the decorator performing pre-processing or post-processing to member functions which will be
   * used
   * by the decorator as a black box, passing through calls which are not overridden.
   */
  std::unique_ptr<FileSampleProvider> fileSampleProvider;
};
} /* namespace datadriven */
} /* namespace sgpp */
