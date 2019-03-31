// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org



#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/PDFFitter.hpp>
#include <sgpp/globaldef.hpp>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <algorithm>
#include <iomanip>

using sgpp::base::Grid;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::HashGridStorage;
using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {


void PDFFitter::fit(Dataset& newDataset, std::unique_ptr<sgpp::base::Grid>& grid2,
        std::vector<size_t> ind, bool val) {
    // Get configurations
    auto datas = newDataset.getData();
    auto& databaseConfig = this->config->getDatabaseConfig();
    auto& gridConfig = this->config->getGridConfig();
    auto& refinementConfig = this->config->getRefinementConfig();
    auto& regularizationConfig = this->config->getRegularizationConfig();
    auto& densityEstimationConfig = this->config->getDensityEstimationConfig();
    // clear model
    reset();
    gridConfig.dim_ = datas.getNcols();
    // TODO(fuchsgruber): Support for geometry aware sparse grids (pass interactions from config?)
    grid = std::unique_ptr<Grid>{buildGrid(gridConfig, &ind)};

    alpha = DataVector{grid->getSize()};
    // Build the offline instance first
    DBMatOffline *offline = nullptr;

    // Intialize database if it is provided
    if (!databaseConfig.filepath.empty()) {
        datadriven::DBMatDatabase database(databaseConfig.filepath);
        // Check if database holds a fitting lhs matrix decomposition
        if (database.hasDataMatrix(gridConfig, refinementConfig, regularizationConfig,
                                   densityEstimationConfig)) {
            std::string offlineFilepath = database.getDataMatrix(gridConfig, refinementConfig,
                    regularizationConfig, densityEstimationConfig);
            offline = DBMatOfflineFactory::buildFromFile(offlineFilepath);
        }
    }
    // Build and decompose offline object if not loaded from database
    if (offline == nullptr) {
        // Build offline object by factory, build matrix and decompose
        offline = DBMatOfflineFactory::buildOfflineObject(gridConfig, refinementConfig,
                regularizationConfig, densityEstimationConfig);
        offline->buildMatrix(grid.get(), regularizationConfig);
        offline->decomposeMatrix(regularizationConfig, densityEstimationConfig);
    }

    online = std::unique_ptr<DBMatOnlineDE>{DBMatOnlineDEFactory::buildDBMatOnlineDE(*offline,
                                            *grid, regularizationConfig.lambda_)};

    online->computeDensityFunction(alpha, datas, *grid,
                                   this->config->getDensityEstimationConfig(), false,
                                   this->config->getCrossvalidationConfig().enable_);


    online->setBeta(this->config->getLearnerConfig().beta);
    online->normalize(alpha, *grid);
        }
    }  // namespace datadriven
}  // namespace sgpp
