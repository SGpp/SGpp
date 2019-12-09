// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/dataSource/NormalizationTransformation.hpp>

#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <valarray>
#include <cmath>
#include <chrono>

namespace sgpp {
namespace datadriven {

NormalizationTransformation::NormalizationTransformation()
    : datasetTransformed(nullptr), datasetInvTransformed(nullptr) {}

void NormalizationTransformation::initialize(Dataset *dataset, DataTransformationConfig config) {
  this->nmConfig = config.normalizationConfig;
  // if MinMax not given by user, p determination
  if (!nmConfig.manualInput) {
    base::DataMatrix& searchData = dataset->getData();
    std::valarray<double> mean(dataset->getDimension());
  // initialize vector
    for (unsigned int d = 0; d < dataset->getDimension(); d++) {
      nmConfig.minmaxInput.push_back({std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()});
    }
    // determine minmax and mean
    // full search
    if (nmConfig.searchInstances == 0) {
        for (size_t i = 0; i < dataset->getNumberInstances() * dataset->getDimension();
             i += dataset->getDimension()) {
               for (unsigned int d = 0; d < dataset->getDimension(); d++) {
                 if (searchData[i+d] < nmConfig.minmaxInput.at(d).at(0)) {
                   nmConfig.minmaxInput.at(d).at(0) = searchData[i+d];
                 }
                 if (searchData[i+d] > nmConfig.minmaxInput.at(d).at(1)) {
                   nmConfig.minmaxInput.at(d).at(1) = searchData[i+d];
                 }
            }
        }
    // partly search with random instances
    } else {
      std::mt19937 generator;
      std::uniform_real_distribution<double> distr(
        0, static_cast<double>(dataset->getNumberInstances() - 1));
      DataVector currSample(dataset->getDimension());
      for (unsigned int i = 0; i < nmConfig.searchInstances; i++) {
        dataset->getData().getRow(static_cast<size_t>(distr(generator)), currSample);
        searchData.setRow(i, currSample);
        }
      for (size_t i = 0; i < nmConfig.searchInstances * dataset->getDimension();
           i += dataset->getDimension()) {
             for (unsigned int d = 0; d < dataset->getDimension(); d++) {
               mean[d] +=searchData[i+d];
               if (searchData[i+d] < nmConfig.minmaxInput.at(d).at(0)) {
                 nmConfig.minmaxInput.at(d).at(0) = searchData[i+d];
               }
               if (searchData[i+d] > nmConfig.minmaxInput.at(d).at(1)) {
                 nmConfig.minmaxInput.at(d).at(1) = searchData[i+d];
               }
             }
      }
    }
    // determine standard deviation or deviation heuristic
    if ((nmConfig.searchInstances != 0) &&
        (nmConfig.searchInstances < dataset->getNumberInstances())) {
          mean = mean / static_cast<double>(nmConfig.searchInstances);
          std::valarray<double> stdDeviation(dataset->getDimension());
          if (nmConfig.stdDeviationHeuristic == false) {
            std::vector<double>stdDeviationTmp(dataset->getDimension());
            for (size_t i = 0; i < nmConfig.searchInstances * dataset->getDimension();
                 i += dataset->getDimension()) {
                   for (unsigned int d = 0; d < dataset->getDimension(); d++) {
                     stdDeviationTmp[d] += ((searchData[i+d] - mean[d]) *
                     (searchData[i+d] - mean[d]));
                   }
            }
            for (unsigned int d = 0; d < dataset->getDimension(); d++) {
              stdDeviation[d] = stdDeviationTmp[d];
            }
            stdDeviation = stdDeviation / double
            (nmConfig.searchInstances - 1);
            stdDeviation = sqrt(stdDeviation);
            } else {
              for (unsigned int d = 0; d < dataset->getDimension(); d++) {
                stdDeviation[d] = (nmConfig.minmaxInput.at(d).at(1) -
                nmConfig.minmaxInput.at(d).at(0)) / nmConfig.deviationHeuristic;
              }
            }
          for (unsigned int d = 0; d < dataset->getDimension(); d++) {
            nmConfig.minmaxInput.at(d).at(0) = nmConfig.minmaxInput.at(d).at(0) -
            nmConfig.minmaxStdDeviation * stdDeviation[d];
            nmConfig.minmaxInput.at(d).at(1) = nmConfig.minmaxInput.at(d).at(1) +
            nmConfig.minmaxStdDeviation * stdDeviation[d];
     }
    }
    std::cout << "Normalization transformation initialized" << std::endl;
  }
}

Dataset *NormalizationTransformation::doTransformation(Dataset *dataset) {
  Dataset datasetTmp = Dataset{dataset->getNumberInstances(), dataset->getDimension()};
  base::DataMatrix& data = dataset->getData();
  base::DataMatrix& dataTransformed = datasetTmp.getData();
  base::DataVector targetsTransformed = base::DataVector();
  base::DataVector& targets = dataset->getTargets();
  base::DataVector tmp;
  size_t outOfBoundC = 0;
  // MinMaxNormalization
  if (nmConfig.method == "minmax") {
    for (size_t i = 0; i < (dataset->getDimension() * dataset->getNumberInstances());
         i += dataset->getDimension()) {
           tmp.resizeZero(0);
           for (unsigned int d = 0; d < dataset->getDimension(); d++) {
             if ((data[i+d] < nmConfig.minmaxInput.at(d).at(0)) ||
                 (data[i+d] > nmConfig.minmaxInput.at(d).at(1))) {
                   outOfBoundC += 1;
                   if (nmConfig.outOfBound == "abort") {
                     std::cout << "Dimension: " << d << std::endl;
                     throw sgpp::base::data_exception("Min or max was not determined correctly");
                   } else if (nmConfig.outOfBound == "discard") {
                     break;
                   } else if (nmConfig.outOfBound == "setToBorder") {
                    if (data[i+d] < nmConfig.minmaxInput.at(d).at(0)) {
                      tmp.append(0);
                    } else {
                      tmp.append(1);
                      }
                   }
             } else {
               tmp.append((data[i+d]  - nmConfig.minmaxInput.at(d).at(0)) /
               (nmConfig.minmaxInput.at(d).at(1) - nmConfig.minmaxInput.at(d).at(0)));
             }
           }
           if (tmp.size() == dataset->getDimension()) {
             dataTransformed.setRow((i - outOfBoundC) / dataset->getDimension(), tmp);
             targetsTransformed.append(targets[i / dataset->getDimension()]);
           }
    }
    std::cout << outOfBoundC << " instances were out of bound." << std::endl;
  }
  if (nmConfig.outOfBound == "discard") {
    dataTransformed.resize(dataset->getNumberInstances() - outOfBoundC);
    datasetTransformed = new Dataset{dataset->getNumberInstances()-outOfBoundC,
                                     dataset->getDimension()};
  } else {
    datasetTransformed = new Dataset{dataset->getNumberInstances(), dataset->getDimension()};
  }
  datasetTransformed->getTargets() = targetsTransformed;
  datasetTransformed->getData() = dataTransformed;
  return datasetTransformed;
}

Dataset *NormalizationTransformation::doInverseTransformation(Dataset *dataset) {
  std::cout << "Performing Denormalization transformation" << std::endl;
  datasetInvTransformed = new Dataset{dataset->getNumberInstances(), dataset->getDimension()};
  base::DataMatrix& data = dataset->getData();
  base::DataMatrix& datatransformed = datasetInvTransformed->getData();
  datasetInvTransformed->getTargets() = dataset->getTargets();

  if (nmConfig.method == "minmax") {
      for (size_t i = 0; i < (dataset->getDimension() * dataset->getNumberInstances());
           i += dataset->getDimension()) {
             for (unsigned int d = 0; d < dataset->getDimension(); d++) {
                  datatransformed[i+d] =
                  data[i+d] * (nmConfig.minmaxInput.at(d).at(1) - nmConfig.minmaxInput.at(d).at(0))
                  + nmConfig.minmaxInput.at(d).at(0);
             }
      }
  }
  return datasetInvTransformed;
}
} /* namespace datadriven */
} /* namespace sgpp */
