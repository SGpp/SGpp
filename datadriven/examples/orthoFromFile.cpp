/*
 * orthoFromFile.cpp
 *
 *  Created on: May 22, 2018
 *      Author: dominik
 */





#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineEigen.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineLU.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>
#include <sgpp/globaldef.hpp>

#include <string>
#include <vector>

int main() {
  sgpp::base::RegularGridConfiguration gridConfig;
    gridConfig.dim_ = 2;
    gridConfig.level_ = 3;
    gridConfig.type_ = sgpp::base::GridType::Linear;

    sgpp::base::AdpativityConfiguration adaptivityConfig;

    sgpp::datadriven::RegularizationConfiguration regularizationConfig;
    regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
    regularizationConfig.lambda_ = 0.1;

    sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
    densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::Eigen;

    sgpp::datadriven::GridFactory gridFactory;
    std::unique_ptr<sgpp::base::Grid> grid = std::unique_ptr<sgpp::base::Grid>{
      gridFactory.createGrid(gridConfig, std::vector<std::vector <size_t>>())
    };

    auto offline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
        sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(gridConfig,
                                                                  adaptivityConfig,
                                                                  regularizationConfig,
                                                                  densityEstimationConfig)};
    offline->buildMatrix(&(*grid), regularizationConfig);
    offline->decomposeMatrix(regularizationConfig, densityEstimationConfig);
    std::string filename = "test.dbmat";
    offline->store(filename);
    std::cout << "stored" << filename << std::endl;
    auto newOffline = std::unique_ptr<sgpp::datadriven::DBMatOffline>{
        sgpp::datadriven::DBMatOfflineFactory::buildFromFile(filename)};

    std::cout << "read" << filename << std::endl;
    std::remove(filename.c_str());

    /**
     * Check matrices
     */
    // lhsMatrix of parent object
    auto& oldMatrix = offline->getDecomposedMatrix();

    std::cout << "Got old matrix decom" << std::endl;
    auto& newMatrix = newOffline->getDecomposedMatrix();

    std::cout << "Got decompositions" << std::endl;
}
