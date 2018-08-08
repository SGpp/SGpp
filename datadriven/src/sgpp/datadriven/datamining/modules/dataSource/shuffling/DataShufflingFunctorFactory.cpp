/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DataShufflingFunctorFactory.cpp
 *
 *  Created on: Jul 23, 2018
 *      Author: dominik
 */

#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorFactory.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorRandom.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorSequential.hpp>

namespace sgpp {
namespace datadriven {
DataShufflingFunctor* DataShufflingFunctorFactory::buildDataShufflingFunctor(
    const DataSourceConfig& config) const {
  switch (config.shuffling) {
    case DataSourceShufflingType::random : {
      return new DataShufflingFunctorRandom(config.randomSeed);
    }
    case DataSourceShufflingType::sequential : {
      return new DataShufflingFunctorSequential();
    }
    default:
      return nullptr;
  }
}
} /* namespace datadriven */
} /* namespace sgpp */




