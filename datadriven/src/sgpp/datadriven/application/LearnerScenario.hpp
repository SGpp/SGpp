// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/base/tools/json/JSON.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

enum class InternalPrecision { Float, Double };

class TestsetConfiguration {
 public:
  bool hasTestDataset;
  std::string alphaReferenceFileName;
  double expectedMSE;
  double expectedLargestDifference;

  TestsetConfiguration()
      : hasTestDataset(false),
        alphaReferenceFileName(""),
        //        valuesFileName(""),
        expectedMSE(0.0),
        expectedLargestDifference(0.0) {}
};

class LearnerScenario : public json::JSON {
 private:
  bool initialized;
  //
  //  // variables for the scenario
  //  double lambda;
  //  std::string datasetFileName;
  //  base::RegularGridConfiguration gridConfig;
  //  solver::SLESolverConfiguration solverConfigRefine;
  //  solver::SLESolverConfiguration solverConfigFinal;
  //  base::AdaptivityConfiguration adaptivityConfig;
  //  datadriven::TestsetConfiguration testsetConfig;

 public:
  LearnerScenario();

  explicit LearnerScenario(std::string scenarioFileName);

  LearnerScenario(std::string datasetFileName, double lambda, InternalPrecision internalPrecision,
                  base::RegularGridConfiguration gridConfig,
                  solver::SLESolverConfiguration SLESolverConfigRefine,
                  solver::SLESolverConfiguration SLESolverConfigFinal,
                  base::AdaptivityConfiguration adaptivityConfig);

  LearnerScenario(std::string datasetFileName, double lambda, InternalPrecision internalPrecision,
                  base::RegularGridConfiguration gridConfig,
                  solver::SLESolverConfiguration SLESolverConfigRefine,
                  solver::SLESolverConfiguration SLESolverConfigFinal,
                  base::AdaptivityConfiguration adaptivityConfig,
                  datadriven::TestsetConfiguration testsetConfiguration);

  bool isInitialized() const;

  void setDatasetFileName(std::string datasetFileName);

  std::string getDatasetFileName();

  void setLambda(double lambda);

  double getLambda();

  void setInternalPrecision(InternalPrecision internalPrecision);

  InternalPrecision getInternalPrecision();

  void setGridConfig(base::RegularGridConfiguration& gridConfig);

  base::RegularGridConfiguration getGridConfig();

  void setSolverConfigurationRefine(solver::SLESolverConfiguration& solverConfigRefine);

  solver::SLESolverConfiguration getSolverConfigurationRefine();

  void setSolverConfigurationFinal(solver::SLESolverConfiguration& solverConfigFinal);

  solver::SLESolverConfiguration getSolverConfigurationFinal();

  void setAdaptivityConfiguration(base::AdaptivityConfiguration& adaptivityConfig);

  base::AdaptivityConfiguration getAdaptivityConfiguration();

  bool hasTestsetConfiguration();

  void setTestsetConfiguration(datadriven::TestsetConfiguration& testsetConfig);

  datadriven::TestsetConfiguration getTestsetConfiguration();

  //  void writeToFile(std::string fileName);
  //
  //  void readFromFile(std::string fileName);

 private:
  template <class T>
  T fromString(const std::string& s);
};
}  // namespace datadriven
}  // namespace sgpp
