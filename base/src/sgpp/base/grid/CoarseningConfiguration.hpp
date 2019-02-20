/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * RefinementConfiguration.hpp
 *
 *  Created on: Aug 4, 2018
 *      Author: dominik
 */

#pragma once

namespace sgpp {
namespace base {
/**
 * Enumeration that defines different types of refinement indicators / functors
 */
enum class CoarseningFunctorType {
  GridPointBased
};

/**
 * Enumeration that defines the different types of refinement monitors (that trigger refinements)
 */
enum class CoarseningMonitorType {
  Error
};
}  // namespace base
}  // namespace sgpp