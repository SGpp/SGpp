#include <iomanip>
#include <iostream>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>

int main(int argc, char** argv) {
  // output path of sereialized test matrix
  std::string outputPath = "tmp";

  // init database with path to hardcoded test database file
  sgpp::datadriven::DBMatDatabase db("database.json");

  // config
  sgpp::base::CombiGridConfiguration gridConfig;

  gridConfig.dim_ = 2;
  gridConfig.levels = std::vector<size_t>{2, 3};

  sgpp::base::AdaptivityConfiguration adpativityConfig;

  adpativityConfig.numRefinements_ = 0;

  sgpp::datadriven::RegularizationConfiguration regConfig;
  regConfig.lambda_ = 0.001;

  sgpp::datadriven::DensityEstimationConfiguration densConfig;
  densConfig.type_ = sgpp::datadriven::DensityEstimationType::Decomposition;
  densConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::OrthoAdapt;

  // create the underlying grid
  std::unique_ptr<sgpp::base::Grid> grid =
      std::unique_ptr<sgpp::base::Grid>{sgpp::base::Grid::createLinearGrid(gridConfig.dim_)};

  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.anisotropicFull(gridConfig.levels);

  sgpp::datadriven::DBMatOfflineOrthoAdapt* mat = (sgpp::datadriven::DBMatOfflineOrthoAdapt*)
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(gridConfig, adpativityConfig,
                                                                regConfig, densConfig);
  // build base matrix
  mat->buildMatrix(grid.get(), regConfig);
  // mat->decomposeMatrix(regConfig, densConfig);

  gridConfig.levels = std::vector<size_t>{3, 2};
  std::unique_ptr<sgpp::base::Grid> grid2 =
      std::unique_ptr<sgpp::base::Grid>{sgpp::base::Grid::createLinearGrid(gridConfig.dim_)};
  sgpp::base::GridGenerator& gridGen2 = grid2->getGenerator();
  gridGen2.anisotropicFull(gridConfig.levels);

  sgpp::datadriven::DBMatOfflineOrthoAdapt* permMat = (sgpp::datadriven::DBMatOfflineOrthoAdapt*)
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(gridConfig, adpativityConfig,
                                                                regConfig, densConfig);
  // build base matrix
  permMat->buildMatrix(grid2.get(), regConfig);
  // permMat->decomposeMatrix(regConfig, densConfig);

  sgpp::base::CombiGridConfiguration baseConfig;
  baseConfig.levels = std::vector<size_t>{2, 3};
  baseConfig.dim_ = 2;

  

  std::cout << "Starting permutation"
            << "\n\n\n";

  std::cout << mat->getLhsMatrix_ONLY_FOR_TESTING().toString() << "\n\n\n";

  std::cout << permMat->getLhsMatrix_ONLY_FOR_TESTING().toString() << "\n\n\n";

  /*std::cout << mat->getTinv().toString() << "\n\n\n";

  std::cout << permMat->getTinv().toString() << "\n\n\n";*/
  mat->permutateLhsMatrix(baseConfig, gridConfig);

  std::cout << mat->getLhsMatrix_ONLY_FOR_TESTING().toString() << "\n\n\n";

  for (size_t i = 0; i < mat->getLhsMatrix_ONLY_FOR_TESTING().getNrows(); i++) {
    for (size_t j = 0; j < mat->getLhsMatrix_ONLY_FOR_TESTING().getNcols(); j++) {
      if (mat->getLhsMatrix_ONLY_FOR_TESTING().get(i, j) -
              permMat->getLhsMatrix_ONLY_FOR_TESTING().get(i, j) >
          0.0001) {
        std::cout << "Unequal pair L: " << i << ", " << j
                  << " mat: " << mat->getLhsMatrix_ONLY_FOR_TESTING().get(i, j)
                  << " perMat: " << permMat->getLhsMatrix_ONLY_FOR_TESTING().get(i, j) << "\n";
      }
    }
  }
  return 0;
}