
/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DataShufflingFunctorSequential.cpp
 *
 *  Created on: Jul 20, 2018
 *      Author: dominik
 */

#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorSequential.hpp>

namespace sgpp {
namespace datadriven {
size_t DataShufflingFunctorSequential::operator()(size_t idx) { return idx; }
} /* namespace datadriven */
} /* namespace sgpp */




