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

#include <memory>

namespace sgpp {
namespace datadriven {

class FileSampleDecorator : public FileSampleProvider {
 public:
  explicit FileSampleDecorator(FileSampleProvider* fileSampleProvider)
      : fileSampleProvider(std::shared_ptr<FileSampleProvider>(fileSampleProvider)) {}
  virtual ~FileSampleDecorator() {}

 protected:
  std::shared_ptr<FileSampleProvider> fileSampleProvider;
};

} /* namespace datadriven */
} /* namespace sgpp */
