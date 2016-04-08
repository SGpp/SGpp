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

#include <sgpp/datadriven/datamining/dataSource/FileSampleProvider.hpp>
#include <memory>

namespace sgpp {
namespace datadriven {

class FileSampleDecorator : public FileSampleProvider {
 public:
  explicit FileSampleDecorator(std::unique_ptr<FileSampleProvider> fileSampleProvider)
      : fileSampleProvider(std::move(fileSampleProvider)) {}
  FileSampleDecorator(FileSampleDecorator&& f)
      : fileSampleProvider(std::move(f.fileSampleProvider)) {}
  virtual ~FileSampleDecorator() {}

 protected:
  std::unique_ptr<FileSampleProvider> fileSampleProvider;
};

} /* namespace datadriven */
} /* namespace sgpp */
