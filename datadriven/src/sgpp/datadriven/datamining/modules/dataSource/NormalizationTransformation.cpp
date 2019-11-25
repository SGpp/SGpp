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

namespace sgpp {
namespace datadriven {

NormalizationTransformation::NormalizationTransformation()
    : datasetTransformed(nullptr), datasetInvTransformed(nullptr) {}


void NormalizationTransformation::initialize(Dataset *dataset, DataTransformationConfig config) {
  this->nmConfig = config.normalizationConfig;
  // if MinMax not given by user, p determination
  if (!nmConfig.manualInput){
	  base::DataMatrix& searchData = dataset->getData();
	  std::valarray<double> mean(dataset->getDimension());
  // initialize vector
    for (unsigned int d = 0; d < dataset->getDimension(); d++){
      nmConfig.minmaxInput.push_back({std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()});
    }
    //determine minmax and mean
    for (unsigned int d = 0; d < dataset->getDimension(); d++){
      for (size_t i = d; i < nmConfig.searchPortion*(dataset->getDimension() *
           dataset->getNumberInstances()); i += dataset->getDimension()){
    	     mean[d] +=searchData[i];
             if (searchData[i] < nmConfig.minmaxInput.at(d).at(0)){
               nmConfig.minmaxInput.at(d).at(0) = searchData[i];
             }
             if (searchData[i] > nmConfig.minmaxInput.at(d).at(1)){
               nmConfig.minmaxInput.at(d).at(1) = searchData[i];
             }
        }
    }
    //determine standard deviation (remove if clause, since only useful if whole data is in dataset
    if (nmConfig.searchPortion != 1){
        mean = mean/(nmConfig.searchPortion * dataset->getNumberInstances());
        std::valarray<double> stdDeviation(dataset->getDimension());

    	for (unsigned int d = 0; d < dataset->getDimension(); d++){
    	  for (size_t i = d; i < nmConfig.searchPortion*(dataset->getDimension() *
               dataset->getNumberInstances()); i += dataset->getDimension()){
    	         stdDeviation[d] += pow(searchData[i], 2);
    	  }
          stdDeviation = stdDeviation / (nmConfig.searchPortion * dataset->getNumberInstances() - 1);
          stdDeviation = sqrt(stdDeviation);
          nmConfig.minmaxInput.at(d).at(0) = nmConfig.minmaxInput.at(d).at(0) -
          nmConfig.minmaxStdDeviation * stdDeviation[d];
          nmConfig.minmaxInput.at(d).at(1) = nmConfig.minmaxInput.at(d).at(1) +
          nmConfig.minmaxStdDeviation * stdDeviation[d];
    	}
    }
  }
  std::cout << "Normalization transformation initialized" << std::endl;
}

Dataset *NormalizationTransformation::doTransformation(Dataset *dataset){
  datasetTransformed = new Dataset{dataset->getNumberInstances(), dataset->getDimension()};
  base::DataMatrix& data = dataset->getData();
  base::DataMatrix& datatransformed = datasetTransformed->getData();
  datasetTransformed->getTargets() = dataset->getTargets();

  //MinMaxNormalization
  if (nmConfig.method == "minmax"){
    for (unsigned int d = 0; d < dataset->getDimension(); d++){
      for (size_t i = d; i < (dataset->getDimension() * dataset->getNumberInstances());
           i += dataset->getDimension()){
             datatransformed[i] = (data[i]  - nmConfig.minmaxInput.at(d).at(0)) /
  	       	 (nmConfig.minmaxInput.at(d).at(1) - nmConfig.minmaxInput.at(d).at(0));
             if ((datatransformed[i] < 0) || (datatransformed[i] > 1)){
               std::cout << "Dimension: " << d << std::endl;
               throw sgpp::base::data_exception("Min or max was not determined correctly");
             }
           }
      }
  }

  return datasetTransformed;
}

Dataset *NormalizationTransformation::doInverseTransformation(Dataset *dataset) {
  std::cout << "Performing Denormalization transformation" << std::endl;
  datasetInvTransformed = new Dataset{dataset->getNumberInstances(), dataset->getDimension()};
  base::DataMatrix& data = dataset->getData();
  base::DataMatrix& datatransformed = datasetTransformed->getData();
  datasetInvTransformed->getTargets() = dataset->getTargets();

  if (nmConfig.method == "minmax"){
      for (unsigned int d = 0; d < dataset->getDimension(); d++){
       for (unsigned int i = d; i < (dataset->getDimension() * dataset->getNumberInstances());
                i += dataset->getDimension()){
                  datatransformed[i] =
                  data[i] * (nmConfig.minmaxInput.at(d).at(1) - nmConfig.minmaxInput.at(d).at(0))
                  + nmConfig.minmaxInput.at(d).at(0);
       }
      }
  }
  return datasetInvTransformed;
}
} /* namespace datadriven */
} /* namespace sgpp */
