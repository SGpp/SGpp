/*
 * test_datadrivenCommon.cpp
 *
 *  Created on: Jul 30, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/globaldef.hpp>

std::string uncompressFile(std::string fileName);

void doRandomRefinements(sgpp::base::AdpativityConfiguration& adaptConfig,
sgpp::base::Grid& grid, sgpp::base::GridGenerator& gridGen);

void doDirectedRefinements(sgpp::base::AdpativityConfiguration& adaptConfig,
sgpp::base::Grid& grid, sgpp::base::GridGenerator& gridGen);
