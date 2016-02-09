// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include "DataMiningConfiguration.hpp"

#include <vector>
#include <map>
#include <string>

namespace SGPP {
namespace datadriven {

class DataMiningConfigurationLeastSquares: public DataMiningConfiguration {
public:
	DataMiningConfigurationLeastSquares();

	explicit DataMiningConfigurationLeastSquares(const std::string& fileName);

//	virtual DataMiningConfiguration* clone() override;
};

}  // namespace base
}  // namespace SGPP

