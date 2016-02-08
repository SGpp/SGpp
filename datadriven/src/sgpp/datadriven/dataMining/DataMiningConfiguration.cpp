// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "DataMiningConfiguration.hpp"

#include <string>

namespace SGPP {
namespace datadriven {

DataMiningConfiguration::DataMiningConfiguration(): json::JSON() {
}

DataMiningConfiguration::DataMiningConfiguration(const std::string& fileName):
  json::JSON(fileName) {
}

DataMiningConfiguration* DataMiningConfiguration::clone() {
	DataMiningConfiguration* clone = new DataMiningConfiguration(*this);
  return clone;
}

}  // namespace base
}  // namespace SGPP

