// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <boost/program_options.hpp>

#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "sgpp/base/datatypes/DataVector.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/grid/GridStorage.hpp"
#include "sgpp/base/grid/generation/GridGenerator.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OpFactory.hpp"
#include "sgpp/datadriven/operation/hash/OperationDensityMultiplicationAVX/OperationDensityMultiplicationAVX.hpp"
#include "sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp"
#include "sgpp/datadriven/operation/hash/OperationPruneGraphOCL/OpFactory.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/globaldef.hpp"
#include "sgpp/solver/sle/ConjugateGradients.hpp"

int main(int argc, char **argv) {
  std::string datasetFileName;
  size_t eval_grid_level;
  size_t level;
  double lambda;

  boost::program_options::options_description description("Allowed options");

  description.add_options()("help", "display help")(
      "datasetFileName", boost::program_options::value<std::string>(&datasetFileName),
      "training data set as an arff file")(
      "eval_grid_level", boost::program_options::value<size_t>(&eval_grid_level)->default_value(2),
      "level for the evaluation of the sparse grid density function (for picture creation)")(
      "level", boost::program_options::value<size_t>(&level)->default_value(4),
      "level of the sparse grid used for density estimation")(
      "lambda", boost::program_options::value<double>(&lambda)->default_value(0.000001),
      "regularization for density estimation");

  boost::program_options::variables_map variables_map;

  boost::program_options::parsed_options options = parse_command_line(argc, argv, description);
  boost::program_options::store(options, variables_map);
  boost::program_options::notify(variables_map);

  if (variables_map.count("help")) {
    std::cout << description << std::endl;
    return 0;
  }

  if (variables_map.count("datasetFileName") == 0) {
    std::cerr << "option \"datasetFileName\" not specified" << std::endl;
    return 1;
  } else {
    std::cout << "datasetFileName: " << datasetFileName << std::endl;
  }

  if (variables_map.count("eval_grid_level") == 0) {
    std::cerr << "option \"eval_grid_level\" not specified" << std::endl;
    return 1;
  } else {
    std::cout << "eval_grid_level: " << eval_grid_level << std::endl;
  }

  if (variables_map.count("level") == 0) {
    std::cerr << "option \"level\" not specified" << std::endl;
    return 1;
  } else {
    std::cout << "level: " << level << std::endl;
  }

  if (variables_map.count("lambda") == 0) {
    std::cerr << "option \"lambda\" not specified" << std::endl;
    return 1;
  } else {
    std::cout << "lambda: " << lambda << std::endl;
  }

  sgpp::datadriven::Dataset dataset =
    sgpp::datadriven::ARFFTools::readARFFFromFile(datasetFileName);
  size_t dimension = dataset.getDimension();
  std::cout << "dimension: " << dimension << std::endl;
  sgpp::base::DataMatrix &trainingData = dataset.getData();

  // Create Grid
  sgpp::base::Grid *grid = sgpp::base::Grid::createLinearGrid(dimension);
  sgpp::base::GridGenerator &gridGen = grid->getGenerator();
  gridGen.regular(level);
  size_t gridsize = grid->getStorage().getSize();
  std::cout << "Grid created! Number of grid points:     " << gridsize << std::endl;

  sgpp::base::DataVector alpha(gridsize);
  sgpp::base::DataVector result(gridsize);
  alpha.setAll(1.0);

  sgpp::solver::ConjugateGradients *solver = new sgpp::solver::ConjugateGradients(1000, 0.0001);
  sgpp::datadriven::DensityOCLMultiPlatform::OperationDensity *operation_mult =
      sgpp::datadriven::createDensityOCLMultiPlatformConfigured(*grid, dimension, lambda,
                                                                "MyOCLConf.cfg");

  std::cout << "Creating rhs" << std::endl;
  sgpp::base::DataVector b(gridsize);
  operation_mult->generateb(trainingData, b);
  // for (size_t i = 0; i < 300; i++) std::cout << b[i] << " ";
  // std::cout << std::endl;
  std::ofstream out_rhs("rhs_erg_dim2_depth11.txt");
  out_rhs.precision(17);
  for (size_t i = 0; i < gridsize; ++i) {
    out_rhs << b[i] << " ";
  }
  out_rhs.close();

  std::cout << "Solving density SLE" << std::endl;
  solver->solve(*operation_mult, alpha, b, false, true);

  std::cout << "Creating evaluation grid" << std::endl;
  double h = 1.0 / std::pow(2.0, eval_grid_level);
  size_t dim_grid_points = 1 << eval_grid_level;
  sgpp::base::DataMatrix evaluationPoints(dim_grid_points * dim_grid_points, 2);
  size_t linearIndex = 0;
  for (size_t i = 0; i < dim_grid_points; i++) {
    for (size_t j = 0; j < dim_grid_points; j++) {
      double x = static_cast<double>(i) * h;
      double y = static_cast<double>(j) * h;
      evaluationPoints(linearIndex, 0) = x;
      evaluationPoints(linearIndex, 1) = y;
      linearIndex += 1;
    }
  }

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
      sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT);

  std::cout << "Creating multieval operation" << std::endl;
  std::unique_ptr<sgpp::base::OperationMultipleEval> eval =
      std::unique_ptr<sgpp::base::OperationMultipleEval>(
          sgpp::op_factory::createOperationMultipleEval(*grid, evaluationPoints, configuration));

  sgpp::base::DataVector results(evaluationPoints.getNrows());
  std::cout << "Evaluating at evaluation grid points" << std::endl;
  eval->mult(alpha, results);

  std::cout << results.toString();
  std::cout << std::endl << "all done!" << std::endl;
}
