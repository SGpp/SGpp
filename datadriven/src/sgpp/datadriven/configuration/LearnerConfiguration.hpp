/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * LearnerConfiguration.hpp
 *
 *  Created on: Jun 19, 2018
 *      Author: dominik
 */

#pragma once

namespace sgpp {
namespace datadriven {

  /**
   * Structure that contains information about the learners behaviour
   */
  struct LearnerConfiguration {
    /**
    * Weigting factor for older batches
    * TODO(fuchsgruber): This is not yet part of DBMatOnlineDE
    */
    double beta = 1.0;
  };
}
}
