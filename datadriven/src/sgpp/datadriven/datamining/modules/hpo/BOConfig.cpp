//
// Created by Eric on 20.03.2018.
//

#include <vector>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <cmath>
#include <iostream>
#include "BOConfig.hpp"
// #include <random>

namespace sgpp {
namespace datadriven {

BOConfig::BOConfig(std::vector<int>* discOptions, std::vector<int>* catOptions, size_t nCont)
: discOptions(discOptions), catOptions(catOptions){
  cont = base::DataVector(nCont, 0);
  disc = std::vector<int>(discOptions->size(), 0);
  cat = std::vector<int>(catOptions->size(), 0);
}

bool BOConfig::nextDisc() {
  for (size_t i = 0; i < disc.size(); ++i) {
    disc[i]++;
    if(disc[i]<discOptions->at(i)){
      return true;
    }
    disc[i] = 0;
  }
  for (size_t i = 0; i < cat.size(); ++i) {
    cat[i]++;
    if(cat[i]<catOptions->at(i)){
      return true;
    }
    cat[i] = 0;
  }
  return false;
}

void BOConfig::calcDiscDistance(BOConfig &other) {
  discDistance = 0;
  for (size_t i = 0; i < disc.size(); ++i) {
    discDistance += std::pow((disc[i] - other.disc[i])/(discOptions->at(i)-1.0),2);
  }
  for (size_t i = 0; i < cat.size(); ++i) {
    discDistance += (cat[i] != other.cat[i]);
  }
}

double BOConfig::getTotalDistance(const base::DataVector &input) {
  double tmp = discDistance;
  for (size_t i = 0; i < cont.size(); ++i) {
    tmp += std::pow((cont[i] - input[i]),2);
  }
  return tmp;
}

double BOConfig::getDistance(BOConfig &other) {
  calcDiscDistance(other);
  return getTotalDistance(other.cont);
}

size_t BOConfig::getContSize() {
  return cont.size();
}

void BOConfig::setCont(const base::DataVector &input) {
  cont = base::DataVector(input);
}

double BOConfig::getCont(int idx) {
  return cont[idx];
}

int BOConfig::getDisc(int idx) {
  return disc[idx];
}

int BOConfig::getCat(int idx) {
  return cat[idx];
}

void BOConfig::setScore(double input) {
  score = input;
}

double BOConfig::getScore() {
  return score;
}

void BOConfig::randomize(std::mt19937& generator) {
  for (size_t i = 0; i < disc.size(); ++i) {
    std::uniform_int_distribution<int> distribution(0, discOptions->at(i)-1);
    disc[i] = distribution(generator);
  }
  for (size_t i = 0; i < cat.size(); ++i) {
    std::uniform_int_distribution<int> distribution(0, catOptions->at(i)-1);
    cat[i] = distribution(generator);
  }
  for (size_t i = 0; i < cont.size(); ++i) {
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    cont[i] = distribution(generator);
  }
}


}
}