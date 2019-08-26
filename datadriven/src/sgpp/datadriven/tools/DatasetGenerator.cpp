// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/tools/DatasetGenerator.hpp>

#include <sgpp/globaldef.hpp>

#include <cstdlib>
#include <cmath>

namespace sgpp {
namespace datadriven {

double DatasetGenerator::uniform(double a, double b) {
  return a + static_cast<double>(rand()) / RAND_MAX * (b - a);
}
DatasetGenerator::~DatasetGenerator() {}

// generate normally distributed random variable
// see http://de.wikipedia.org/wiki/Box-Muller-Methode
/**
 * @brief generate normally distributed random variable
 *        see http://de.wikipedia.org/wiki/Box-Muller-Methode
 * uses rand(), make sure to call srand yourself
 * @param mean mean of distribution
 * @param stddev std deviation of distribution
 * @return normally distributed random var
 */
double DatasetGenerator::normal(double mean, double stddev) {
  double u1 = static_cast<double>(rand()) / RAND_MAX;
  double u2 = static_cast<double>(rand()) / RAND_MAX;
  double stdnormal = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
  return mean + stdnormal * stddev;
}

void Friedman1Generator::createData(size_t offset, size_t size,
                                    base::DataMatrix& trData, base::DataVector& classes) {
  trData.resize(size);
  classes.resize(size);
  srand(static_cast<unsigned int>(offset));

  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < 10; j++) {
      trData.set(i, j, uniform(0.0, 1.0));
    }

    double cls = 10.0 * sin(M_PI * trData.get(i, 0) * trData.get(i, 1))
                  + 20.0 * (trData.get(i, 2) - 0.5) * (trData.get(i, 2) - 0.5)
                  + 10.0 * trData.get(i, 3)
                  + 5.0 * trData.get(i, 4) + normal(0.0, 1.0);
    classes.set(i, cls);
  }
}

size_t Friedman1Generator::getDims() {
  return 10;
}



void Friedman2Generator::createData(size_t offset, size_t size,
                                    base::DataMatrix& trData, base::DataVector& classes) {
  trData.resize(size);
  classes.resize(size);
  srand(static_cast<unsigned int>(offset));

  for (size_t i = 0; i < size; i++) {
    trData.set(i, 0, uniform(0.0, 100.0));
    trData.set(i, 1, uniform(40.0 * M_PI, 560.0 * M_PI));
    trData.set(i, 2, uniform(0.0, 1.0));
    trData.set(i, 3, uniform(1.0, 11.0));
    double tmp = trData.get(i, 1) * trData.get(i, 2) - 1.0 / (trData.get(i, 1) * trData.get(i, 3));
    double cls = sqrt(trData.get(i, 0) * trData.get(i, 0) + tmp * tmp) + normal(0.0, 125.0);
    classes.set(i, cls);
  }

  for (size_t i = 0; i < 4; i++) {
    trData.normalizeDimension(i);
  }
}

size_t Friedman2Generator::getDims() {
  return 4;
}


void Friedman3Generator::createData(size_t offset, size_t size,
                                    base::DataMatrix& trData, base::DataVector& classes) {
  trData.resize(size);
  classes.resize(size);
  srand(static_cast<unsigned int>(offset));

  for (size_t i = 0; i < size; i++) {
    trData.set(i, 0, uniform(0.0, 100.0));
    trData.set(i, 1, uniform(40.0 * M_PI, 560.0 * M_PI));
    trData.set(i, 2, uniform(0.0, 1.0));
    trData.set(i, 3, uniform(1.0, 11.0));
    double tmp = trData.get(i, 1) * trData.get(i, 2) - 1.0 / (trData.get(i,
                  1) * trData.get(i, 3));
    double cls = atan(tmp / trData.get(i, 0)) + normal(0.0, 0.1);
    classes.set(i, cls);
  }

  for (size_t i = 0; i < 4; i++) {
    trData.normalizeDimension(i);
  }
}

size_t Friedman3Generator::getDims() {
  return 4;
}

}  // namespace datadriven
}  // namespace sgpp

