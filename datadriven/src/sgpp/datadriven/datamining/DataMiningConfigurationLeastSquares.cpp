// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "DataMiningConfigurationLeastSquares.hpp"

#include <string>

namespace SGPP {
namespace datadriven {

DataMiningConfigurationLeastSquares::DataMiningConfigurationLeastSquares() :
		DataMiningConfiguration() {
}

DataMiningConfigurationLeastSquares::DataMiningConfigurationLeastSquares(
		const std::string& fileName) :
		DataMiningConfiguration(fileName) {
}

//DataMiningConfiguration* DataMiningConfiguration::clone() {
	//TODO: implement
//	DataMiningConfiguration* clone = new DataMiningConfigurationLeastSquares(*this);
//	return clone;
//}

}  // namespace base
}  // namespace SGPP

