/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FileSampleProviderBuilder.cpp
 *
 *  Created on: 16.05.2016
 *      Author: Michael Lettrich
 */

#include "DataSourceBuilder.hpp"

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/ArffFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/GzipFileSampleDecorator.hpp>

#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {

DataSourceBuilder::DataSourceBuilder() : batchSize(0), numBatches(0) {}

DataSourceBuilder::~DataSourceBuilder() {}

} /* namespace datadriven */
} /* namespace sgpp */
