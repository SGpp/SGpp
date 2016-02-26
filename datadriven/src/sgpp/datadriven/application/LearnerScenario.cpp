// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <string>
#include <vector>

#include "sgpp/datadriven/application/LearnerScenario.hpp"
#include "sgpp/base/exception/not_implemented_exception.hpp"

namespace SGPP {
namespace datadriven {

LearnerScenario::LearnerScenario() : isInitialized(false) {}

LearnerScenario::LearnerScenario(std::string scenarioFileName)
    : json::JSON(scenarioFileName), isInitialized(true) {
  //  this->readFromFile(scenarioFileName);
}

LearnerScenario::LearnerScenario(std::string datasetFileName, double lambda,
                                 InternalPrecision internalPrecision,
                                 base::RegularGridConfiguration gridConfig,
                                 solver::SLESolverConfiguration SLESolverConfigRefine,
                                 solver::SLESolverConfiguration SLESolverConfigFinal,
                                 base::AdpativityConfiguration adaptConfig)
    : isInitialized(true) {
  this->setDatasetFileName(datasetFileName);
  this->setLambda(lambda);
  this->setInternalPrecision(internalPrecision);
  this->setGridConfig(gridConfig);
  this->setSolverConfigurationRefine(SLESolverConfigRefine);
  this->setSolverConfigurationFinal(SLESolverConfigFinal);
  this->setAdaptivityConfiguration(adaptConfig);
  (*this).addDictAttr("testset").addIDAttr("hasTestDataset", false);
}

LearnerScenario::LearnerScenario(std::string datasetFileName, double lambda,
                                 InternalPrecision internalPrecision,
                                 base::RegularGridConfiguration gridConfig,
                                 solver::SLESolverConfiguration SLESolverConfigRefine,
                                 solver::SLESolverConfiguration SLESolverConfigFinal,
                                 base::AdpativityConfiguration adaptConfig,
                                 datadriven::TestsetConfiguration testsetConfig)
    : isInitialized(true) {
  this->setDatasetFileName(datasetFileName);
  this->setLambda(lambda);
  this->setInternalPrecision(internalPrecision);
  this->setGridConfig(gridConfig);
  this->setSolverConfigurationRefine(SLESolverConfigRefine);
  this->setSolverConfigurationFinal(SLESolverConfigFinal);
  this->setAdaptivityConfiguration(adaptConfig);
  this->setTestsetConfiguration(testsetConfig);
}

void LearnerScenario::setDatasetFileName(std::string datasetFileName) {
  (*this).replaceTextAttr("datasetFileName", datasetFileName);
}

std::string LearnerScenario::getDatasetFileName() { return (*this)["datasetFileName"].get(); }

void LearnerScenario::setLambda(double lambda) { (*this).replaceIDAttr("lambda", lambda); }

double LearnerScenario::getLambda() { return (*this)["lambda"].getDouble(); }

void LearnerScenario::setInternalPrecision(InternalPrecision internalPrecision) {
  if (internalPrecision == InternalPrecision::Float) {
    (*this).replaceIDAttr("INTERNAL_PRECISION", "float");
  } else {
    (*this).replaceIDAttr("INTERNAL_PRECISION", "double");
  }
}

InternalPrecision LearnerScenario::getInternalPrecision() {
  if ((*this)["INTERNAL_PRECISION"].get().compare("float") == 0) {
    return InternalPrecision::Float;
  } else if ((*this)["INTERNAL_PRECISION"].get().compare("double") == 0) {
    return InternalPrecision::Double;
  } else {
    throw base::not_implemented_exception(
        "error: learner does only support \"float\" and \"double\"");
  }
}

void LearnerScenario::setGridConfig(base::RegularGridConfiguration& gridConfig) {
  (*this).replaceDictAttr("grid");
  (*this)["grid"].replaceIDAttr("boundaryLevel", gridConfig.boundaryLevel_);
  (*this)["grid"].replaceIDAttr("dim", gridConfig.dim_);
  (*this)["grid"].replaceIDAttr("level", static_cast<uint64_t>(gridConfig.level_));
  (*this)["grid"].replaceIDAttr("maxDegree", gridConfig.maxDegree_);

  if (gridConfig.type_ == base::GridType::Linear) {
    (*this)["grid"].replaceTextAttr("type", "Linear");
  } else if (gridConfig.type_ == base::GridType::ModLinear) {
    (*this)["grid"].replaceTextAttr("type", "ModLinear");
  } else {
    throw base::not_implemented_exception(
        "error: learner does not support the specified grid type");
  }
}

base::RegularGridConfiguration LearnerScenario::getGridConfig() {
  base::RegularGridConfiguration gridConfig;
  gridConfig.boundaryLevel_ = (*this)["grid"]["boundaryLevel"].getUInt();
  gridConfig.dim_ = (*this)["grid"]["dim"].getUInt();
  gridConfig.level_ = static_cast<int>((*this)["grid"]["level"].getInt());
  gridConfig.maxDegree_ = (*this)["grid"]["maxDegree"].getUInt();

  std::string typeString = (*this)["grid"]["type"].get();
  if (typeString.compare("Linear") == 0) {
    gridConfig.type_ = base::GridType::Linear;
  } else if (typeString.compare("ModLinear") == 0) {
    gridConfig.type_ = base::GridType::ModLinear;
  } else {
    throw base::not_implemented_exception(
        "error: learner does not support the specified grid type");
  }
  return gridConfig;
}

void LearnerScenario::setSolverConfigurationRefine(
    solver::SLESolverConfiguration& solverConfigRefine) {
  (*this).replaceDictAttr("solverRefine");
  (*this)["solverRefine"].replaceIDAttr("eps", static_cast<double>(solverConfigRefine.eps_));
  (*this)["solverRefine"].replaceIDAttr("maxIterations", solverConfigRefine.maxIterations_);
  (*this)["solverRefine"].replaceIDAttr("threshold", solverConfigRefine.threshold_);
  if (solverConfigRefine.type_ == solver::SLESolverType::CG) {
    (*this)["solverRefine"].replaceIDAttr("type", "CG");
  } else if (solverConfigRefine.type_ == solver::SLESolverType::BiCGSTAB) {
    (*this)["solverRefine"].replaceIDAttr("type", "BiCGSTAB");
  } else {
    throw base::not_implemented_exception(
        "error: learner does not support the specified solver type");
  }
}

solver::SLESolverConfiguration LearnerScenario::getSolverConfigurationRefine() {
  solver::SLESolverConfiguration solverConfigFinal;
  solverConfigFinal.eps_ = (*this)["solverRefine"]["eps"].getDouble();
  solverConfigFinal.maxIterations_ = (*this)["solverRefine"]["maxIterations"].getUInt();
  solverConfigFinal.threshold_ = (*this)["solverRefine"]["threshold"].getDouble();
  std::string solverType = (*this)["solverRefine"]["type"].get();
  if (solverType.compare("CG") == 0) {
    solverConfigFinal.type_ = solver::SLESolverType::CG;
  } else if (solverType.compare("BiCGSTAB") == 0) {
    solverConfigFinal.type_ = solver::SLESolverType::BiCGSTAB;
  } else {
    throw base::not_implemented_exception(
        "error: learner does not support the specified solver type");
  }
  return solverConfigFinal;
}

void LearnerScenario::setSolverConfigurationFinal(
    solver::SLESolverConfiguration& solverConfigFinal) {
  (*this).replaceDictAttr("solverFinal");
  (*this)["solverFinal"].replaceIDAttr("eps", static_cast<double>(solverConfigFinal.eps_));
  (*this)["solverFinal"].replaceIDAttr("maxIterations", solverConfigFinal.maxIterations_);
  (*this)["solverFinal"].replaceIDAttr("threshold", solverConfigFinal.threshold_);
  if (solverConfigFinal.type_ == solver::SLESolverType::CG) {
    (*this)["solverFinal"].replaceIDAttr("type", "CG");
  } else if (solverConfigFinal.type_ == solver::SLESolverType::BiCGSTAB) {
    (*this)["solverFinal"].replaceIDAttr("type", "BiCGSTAB");
  } else {
    throw base::not_implemented_exception(
        "error: learner does not support the specified solver type");
  }
}

solver::SLESolverConfiguration LearnerScenario::getSolverConfigurationFinal() {
  solver::SLESolverConfiguration solverConfigFinal;
  solverConfigFinal.eps_ = (*this)["solverFinal"]["eps"].getDouble();
  solverConfigFinal.maxIterations_ = (*this)["solverFinal"]["maxIterations"].getUInt();
  solverConfigFinal.threshold_ = (*this)["solverFinal"]["threshold"].getDouble();
  std::string solverType = (*this)["solverFinal"]["type"].get();
  if (solverType.compare("CG") == 0) {
    solverConfigFinal.type_ = solver::SLESolverType::CG;
  } else if (solverType.compare("BiCGSTAB") == 0) {
    solverConfigFinal.type_ = solver::SLESolverType::BiCGSTAB;
  } else {
    throw base::not_implemented_exception(
        "error: learner does not support the specified solver type");
  }
  return solverConfigFinal;
}

void LearnerScenario::setAdaptivityConfiguration(base::AdpativityConfiguration& adaptConfig) {
  (*this).replaceDictAttr("adaptivity");
  (*this)["adaptivity"].replaceIDAttr("maxLevelType", adaptConfig.maxLevelType_);
  (*this)["adaptivity"].replaceIDAttr("noPoints", adaptConfig.noPoints_);
  (*this)["adaptivity"].replaceIDAttr("numRefinements", adaptConfig.numRefinements_);
  (*this)["adaptivity"].replaceIDAttr("percent", adaptConfig.percent_);
  (*this)["adaptivity"].replaceIDAttr("threshold", adaptConfig.threshold_);
}

base::AdpativityConfiguration LearnerScenario::getAdaptivityConfiguration() {
  base::AdpativityConfiguration adaptConfig;
  adaptConfig.maxLevelType_ = (*this)["adaptivity"]["maxLevelType"].getBool();
  adaptConfig.noPoints_ = (*this)["adaptivity"]["noPoints"].getUInt();
  adaptConfig.numRefinements_ = (*this)["adaptivity"]["numRefinements"].getUInt();
  adaptConfig.percent_ = (*this)["adaptivity"]["percent"].getDouble();
  adaptConfig.threshold_ = (*this)["adaptivity"]["threshold"].getDouble();
  return adaptConfig;
}

template <class T>
T LearnerScenario::fromString(const std::string& s) {
  std::istringstream stream(s);
  T t;
  stream >> t;
  return t;
}

bool LearnerScenario::hasTestsetConfiguration() {
  return (*this)["testset"]["hasTestDataset"].getBool();
}

void LearnerScenario::setTestsetConfiguration(datadriven::TestsetConfiguration& testsetConfig) {
  (*this).replaceDictAttr("testset");
  (*this)["testset"].replaceIDAttr("hasTestDataset", testsetConfig.hasTestDataset);
  (*this)["testset"].replaceTextAttr("testFileName", testsetConfig.alphaReferenceFileName);
  (*this)["testset"].replaceIDAttr("expectedMSE", testsetConfig.expectedMSE);
  (*this)["testset"].replaceIDAttr("expectedLargestDifference",
                                   testsetConfig.expectedLargestDifference);
}

datadriven::TestsetConfiguration LearnerScenario::getTestsetConfiguration() {
  TestsetConfiguration testsetConfig;
  testsetConfig.hasTestDataset = (*this)["testset"]["hasTestDataset"].getBool();
  testsetConfig.alphaReferenceFileName = (*this)["testset"]["testFileName"].get();
  testsetConfig.expectedMSE = (*this)["testset"]["expectedMSE"].getDouble();
  testsetConfig.expectedLargestDifference =
      (*this)["testset"]["expectedLargestDifference"].getDouble();
  return testsetConfig;
}

}  // namespace datadriven
}  // namespace SGPP
