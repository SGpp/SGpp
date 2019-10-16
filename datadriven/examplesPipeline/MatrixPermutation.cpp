
#include <iomanip>
#include <iostream>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>

int main(int argc, char** argv) {
  // output path of sereialized test matrix
  std::string outputPath = "tmp";

  // init database with path to hardcoded test database file
  sgpp::datadriven::DBMatDatabase db("database.json");

  // config
  sgpp::base::CombiGridConfiguration baseGridConfig;

  baseGridConfig.dim_ = 4;
  baseGridConfig.levels = std::vector<size_t>{3, 3, 2, 3};

  sgpp::base::CombiGridConfiguration desiredGridConfig;
  desiredGridConfig.dim_ = 8;
  desiredGridConfig.levels = std::vector<size_t>{1, 1, 2, 1, 1, 3, 3, 3};

  sgpp::base::AdaptivityConfiguration adpativityConfig;

  adpativityConfig.numRefinements_ = 0;

  sgpp::datadriven::RegularizationConfiguration regConfig;
  regConfig.lambda_ = 0;

  sgpp::datadriven::DensityEstimationConfiguration densConfig;
  densConfig.type_ = sgpp::datadriven::DensityEstimationType::Decomposition;
  densConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::OrthoAdapt;

  // create the underlying grid
  std::unique_ptr<sgpp::base::Grid> grid =
      std::unique_ptr<sgpp::base::Grid>{sgpp::base::Grid::createLinearGrid(baseGridConfig.dim_)};

  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
  gridGen.anisotropicFull(baseGridConfig.levels);

  sgpp::datadriven::DBMatOfflineOrthoAdapt* baseMat = (sgpp::datadriven::DBMatOfflineOrthoAdapt*)
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(baseGridConfig, adpativityConfig,
                                                                regConfig, densConfig);

  std::cout << "Building base matrix"
            << "\n";
  // build base matrix
  baseMat->buildMatrix(grid.get(), regConfig);

  std::unique_ptr<sgpp::base::Grid> grid2 =
      std::unique_ptr<sgpp::base::Grid>{sgpp::base::Grid::createLinearGrid(desiredGridConfig.dim_)};
  sgpp::base::GridGenerator& gridGen2 = grid2->getGenerator();
  gridGen2.anisotropicFull(desiredGridConfig.levels);

  sgpp::datadriven::DBMatOfflineOrthoAdapt* permMat = (sgpp::datadriven::DBMatOfflineOrthoAdapt*)
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(desiredGridConfig, adpativityConfig,
                                                                regConfig, densConfig);

  std::cout << "Building perm matrix"
            << "\n";
  // build base matrix
  permMat->buildMatrix(grid2.get(), regConfig);
  // permMat->decomposeMatrix(regConfig, densConfig);

  std::cout << "Decomposing base matrix"
            << "\n";
  baseMat->decomposeMatrix(regConfig, densConfig);

  std::cout << "Decomposing perm matrix"
            << "\n";
  
  permMat->decomposeMatrix(regConfig, densConfig);

  std::cout << "Permutating base matrix"
            << "\n";         
  baseMat->permutateDecomposition(baseGridConfig, desiredGridConfig);



  // Test online fitting
  std::unique_ptr<sgpp::datadriven::DBMatOnlineDE> online1{
      sgpp::datadriven::DBMatOnlineDEFactory::buildDBMatOnlineDE(*baseMat, *grid2,
                                                                 regConfig.lambda_)};

  std::unique_ptr<sgpp::datadriven::DBMatOnlineDE> online2{
      sgpp::datadriven::DBMatOnlineDEFactory::buildDBMatOnlineDE(*permMat, *grid2,
                                                                 regConfig.lambda_)};

  // Generate sample dataset
  sgpp::base::DataMatrix samples(1000, desiredGridConfig.dim_);
  for (int i = 0; i < 1000; i++) {
    sgpp::base::DataVector vec(desiredGridConfig.dim_);
    for (int j = 0; j < desiredGridConfig.dim_; j++) {
      vec.at(j) = (double)std::rand() / RAND_MAX;
    }
    samples.setRow(i, vec);
  }

  sgpp::base::DataVector alpha1(baseMat->getGridSize());
  sgpp::base::DataVector alpha2(permMat->getGridSize());

  online1->computeDensityFunction(alpha1, samples, *grid2, densConfig, false);
  online2->computeDensityFunction(alpha2, samples, *grid2, densConfig, false);

  online1->setBeta(1);
  online2->setBeta(1);

  if (alpha1.getSize() != alpha2.getSize())
    std::cout << "Unequal alpha size."
              << "\n";

  for (int i = 0; i < alpha1.getSize(); i++) {
    if (alpha1[i] - alpha2[i] > 0.001) {
      std::cout << "Uneq alpha: " << alpha1[i] << " | " << alpha2[i] << "\n";
    }
  }

  std::cout << "Starting test"
            << "\n";
  for (int i = 0; i < 10; i++) {
    sgpp::base::DataVector p(desiredGridConfig.dim_);

    for (int j = 0; j < desiredGridConfig.dim_; j++) {
      p.at(j) = (double)std::rand() / RAND_MAX;
    }

    double eval1 = online1->eval(alpha1, p, *grid2);
    double eval2 = online2->eval(alpha2, p, *grid2);

    if (eval1 - eval2 > 0.001) {
      std::cout << "Unequal eval. Point: (" << p.at(0) << "," << p.at(1) << ") Values: " << eval1
                << " | " << eval2 << "\n";
    }
  }

  sgpp::datadriven::DBMatOfflineOrthoAdapt* baseMat1 = (sgpp::datadriven::DBMatOfflineOrthoAdapt*)
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(baseGridConfig, adpativityConfig,
                                                                regConfig, densConfig);

  sgpp::datadriven::DBMatOfflineOrthoAdapt* permMat1 = (sgpp::datadriven::DBMatOfflineOrthoAdapt*)
      sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(desiredGridConfig, adpativityConfig,
                                                                regConfig, densConfig);

  baseMat1->buildMatrix(grid.get(), regConfig);
  permMat1->buildMatrix(grid2.get(), regConfig);

  baseMat1->permutateLhsMatrix(baseGridConfig, desiredGridConfig);

  for (int i = 0; i < permMat1->getGridSize(); i++) {
    for (int j = 0; j < permMat1->getGridSize(); j++) {
      if (baseMat1->getLhsMatrix_ONLY_FOR_TESTING().get(i, j) -
              permMat1->getLhsMatrix_ONLY_FOR_TESTING().get(i, j) >
          0.001) {
        std::cout << "Uneq, base: " << baseMat1->getLhsMatrix_ONLY_FOR_TESTING().get(i, j)
                  << "  perm: " << permMat1->getLhsMatrix_ONLY_FOR_TESTING().get(i, j) << "\n";
      }
    }
  }

  return 0;
}