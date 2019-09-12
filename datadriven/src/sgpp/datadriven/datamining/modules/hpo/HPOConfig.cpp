// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/hpo/HPOConfig.hpp>

#include <vector>


namespace sgpp {
namespace datadriven {

void HPOConfig::setupDefaults() {
  seed = 42;
  nTrainSamples = -1;
  stages = {200, 200, 100};
  constraints = {2, 2};
  lambda = 1;
  nRandom = 10;
}

int64_t HPOConfig::getSeed() const {
  return seed;
}

void HPOConfig::setSeed(int64_t seed) {
  HPOConfig::seed = seed;
}

const std::vector<int64_t> &HPOConfig::getStages() const {
  return stages;
}

void HPOConfig::setStages(const std::vector<int64_t> &stages) {
  HPOConfig::stages = stages;
}

const std::vector<int64_t> &HPOConfig::getConstraints() const {
  return constraints;
}

void HPOConfig::setConstraints(const std::vector<int64_t> &constraints) {
  HPOConfig::constraints = constraints;
}

double HPOConfig::getLambda() const {
  return lambda;
}

void HPOConfig::setLambda(double lambda) {
  HPOConfig::lambda = lambda;
}

int64_t HPOConfig::getNRandom() const {
  return nRandom;
}

void HPOConfig::setNRandom(int64_t nRandom) {
  HPOConfig::nRandom = nRandom;
}

int64_t HPOConfig::getNRuns() const {
  return nRuns;
}

void HPOConfig::setNRuns(int64_t nRuns) {
  HPOConfig::nRuns = nRuns;
}

int64_t HPOConfig::getNTrainSamples() const {
  return nTrainSamples;
}

void HPOConfig::setNTrainSamples(int64_t nTrainSamples) {
  HPOConfig::nTrainSamples = nTrainSamples;
}
} /* namespace datadriven */
} /* namespace sgpp */
