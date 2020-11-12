// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/globaldef.hpp>

#include <string>

std::string uncompressFile(std::string fileName);

void doRandomRefinements(sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                         sgpp::base::Grid& grid, sgpp::base::GridGenerator& gridGen);

void doDirectedRefinements(sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                           sgpp::base::Grid& grid, sgpp::base::GridGenerator& gridGen);
