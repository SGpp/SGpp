// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <string>

#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/globaldef.hpp"

std::string uncompressFile(std::string fileName);

void doRandomRefinements(SGPP::base::AdpativityConfiguration& adaptConfig, SGPP::base::Grid& grid,
                         SGPP::base::GridGenerator& gridGen);

void doDirectedRefinements(SGPP::base::AdpativityConfiguration& adaptConfig, SGPP::base::Grid& grid,
                           SGPP::base::GridGenerator& gridGen);
