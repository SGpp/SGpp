/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * BOConfig.cpp
 *
 *  Created on:	20.03.2018
 *      Author: Eric Koepke
 */
#include <sgpp/datadriven/datamining/modules/hpo/bo/BOConfig.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>

#include <vector>
#include <cmath>
#include <iostream>

namespace sgpp {
namespace datadriven {

BOConfig::BOConfig(std::vector<int> *discOptions, std::vector<int> *catOptions, size_t nCont)
    : discOptions(discOptions), catOptions(catOptions) {
  cont = base::DataVector(nCont, 0);
  disc = std::vector<int>(discOptions->size(), 0);
  cat = std::vector<int>(catOptions->size(), 0);
}

bool BOConfig::nextDisc() {
  for (size_t i = 0; i < disc.size(); ++i) {
    disc[i]++;
    if (disc[i] < discOptions->at(i)) {
      return true;
    }
    disc[i] = 0;
  }
  for (size_t i = 0; i < cat.size(); ++i) {
    cat[i]++;
    if (cat[i] < catOptions->at(i)) {
      return true;
    }
    cat[i] = 0;
  }
  return false;
}

void BOConfig::calcDiscDistance(BOConfig &other, base::DataVector &scales) {
  discDistance = 0;
  size_t k = cont.size();
  for (size_t i = 0; i < disc.size(); ++i) {
    discDistance += std::pow(scales[k] * (disc[i] - other.disc[i]) / (discOptions->at(i) - 1.0), 2);
    k++;
  }
  for (size_t i = 0; i < cat.size(); ++i) {
    discDistance += std::pow(scales[k] * (cat[i] != other.cat[i]), 2);
    k++;
  }
}

double BOConfig::getTotalDistance(const base::DataVector &input, base::DataVector &scales) {
  double tmp = discDistance;
  size_t k = 0;
  for (size_t i = 0; i < cont.size(); ++i) {
    tmp += std::pow(scales[k] * (cont[i] - input[i]), 2);
    k++;
  }
  return tmp;
}

size_t BOConfig::getContSize() {
  return cont.size();
}

size_t BOConfig::getNPar() const {
  return cont.size() + disc.size() + cat.size();
}

void BOConfig::setCont(const base::DataVector &input) {
  cont = base::DataVector(input);
}

double BOConfig::getCont(size_t idx) {
  return cont[idx];
}

int BOConfig::getDisc(size_t idx) {
  return disc[idx];
}

int BOConfig::getCat(size_t idx) {
  return cat[idx];
}

void BOConfig::setScore(double input) {
  score = input;
}

double BOConfig::getScore() {
  return score;
}

void BOConfig::randomize(std::mt19937 &generator) {
  for (size_t i = 0; i < disc.size(); ++i) {
    std::uniform_int_distribution<int> distribution(0, discOptions->at(i) - 1);
    disc[i] = distribution(generator);
  }
  for (size_t i = 0; i < cat.size(); ++i) {
    std::uniform_int_distribution<int> distribution(0, catOptions->at(i) - 1);
    cat[i] = distribution(generator);
  }
  for (size_t i = 0; i < cont.size(); ++i) {
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    cont[i] = distribution(generator);
  }
}

double BOConfig::getScaledDistance(BOConfig &other, const base::DataVector &scales) {
  double tmp = 0;
  size_t k = 0;
  for (size_t i = 0; i < cont.size(); ++i) {
    tmp += std::pow(scales[k] * (cont[i] - other.cont[i]), 2);
    k++;
  }
  for (size_t i = 0; i < disc.size(); ++i) {
    tmp += std::pow(scales[k] * (disc[i] - other.disc[i]) / (discOptions->at(i) - 1.0), 2);
    k++;
  }
  for (size_t i = 0; i < cat.size(); ++i) {
    tmp += std::pow(scales[k] * (cat[i] != other.cat[i]), 2);
    k++;
  }
  return tmp;
}
} /* namespace datadriven */
} /* namespace sgpp */
