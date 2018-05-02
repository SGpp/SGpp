/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * HPOConfig.hpp
 *
 *  Created on: 01.05.2019
 *      Author: Eric Koepke
 */

#pragma once


#include <cstdint>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * Configuration for fitter scenarios using least squares optimization.
 */
class HPOConfig {
 public:
  HPOConfig() = default;


  void setupDefaults();


  int64_t getSeed() const;

  void setSeed(int64_t seed);

  const std::vector<int64_t> &getStages() const;

  void setStages(const std::vector<int64_t> &stages);

  const std::vector<int64_t> &getConstraints() const;

  void setConstraints(const std::vector<int64_t> &constraints);

  double getLambda() const;

  void setLambda(double lambda);

  int64_t getNRandom() const;

  void setNRandom(int64_t nRandom);

  int64_t getNRuns() const;

  void setNRuns(int64_t nRuns);

private:
  int64_t seed;
  std::vector<int64_t> stages;
  std::vector<int64_t> constraints;
  double lambda;
  int64_t nRandom;
  int64_t nRuns;
};

} /* namespace datadriven */
} /* namespace sgpp */
