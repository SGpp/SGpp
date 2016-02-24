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
                                 base::RegularGridConfiguration gridConfig,
                                 solver::SLESolverConfiguration SLESolverConfigRefine,
                                 solver::SLESolverConfiguration SLESolverConfigFinal,
                                 base::AdpativityConfiguration adaptConfig)
    : isInitialized(true) {
  this->setDatasetFileName(datasetFileName);
  this->setLambda(lambda);
  this->setGridConfig(gridConfig);
  this->setSolverConfigurationRefine(SLESolverConfigRefine);
  this->setSolverConfigurationFinal(SLESolverConfigFinal);
  this->setAdaptivityConfiguration(adaptConfig);
  (*this).addDictAttr("testset").addIDAttr("hasTestDataset", false);
}

LearnerScenario::LearnerScenario(std::string datasetFileName, double lambda,
                                 base::RegularGridConfiguration gridConfig,
                                 solver::SLESolverConfiguration SLESolverConfigRefine,
                                 solver::SLESolverConfiguration SLESolverConfigFinal,
                                 base::AdpativityConfiguration adaptConfig,
                                 datadriven::TestsetConfiguration testsetConfig)
    : isInitialized(true) {
  this->setDatasetFileName(datasetFileName);
  this->setLambda(lambda);
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

void LearnerScenario::setGridConfig(base::RegularGridConfiguration& gridConfig) {
  (*this).replaceDictAttr("grid");
  (*this)["grid"].replaceIDAttr("boundaryLevel", gridConfig.boundaryLevel_);
  (*this)["grid"].replaceIDAttr("dim", gridConfig.dim_);
  (*this)["grid"].replaceIDAttr("level", static_cast<uint64_t>(gridConfig.level_));
  (*this)["grid"].replaceIDAttr("maxDegree", gridConfig.maxDegree_);

  if (gridConfig.type_ == base::GridType::Linear) {
    (*this)["grid"].replaceTextAttr("type", "Linear");
  } else if (gridConfig.type_ == base::GridType::Linear) {
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
  }
}

solver::SLESolverConfiguration LearnerScenario::getSolverConfigurationRefine() {
  solver::SLESolverConfiguration solverConfigFinal;
  solverConfigFinal.eps_ = (*this)["solverRefine"]["eps"].getDouble();
  solverConfigFinal.maxIterations_ = (*this)["solverRefine"]["maxIterations"].getUInt();
  solverConfigFinal.threshold_ = (*this)["solverRefine"]["threshold"].getDouble();
  if ((*this)["solverRefine"]["type"].get().compare("CG")) {
    solverConfigFinal.type_ = solver::SLESolverType::CG;
  } else if ((*this)["solverRefine"]["type"].get().compare("BiCGSTAB")) {
    solverConfigFinal.type_ = solver::SLESolverType::BiCGSTAB;
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
  }
}

solver::SLESolverConfiguration LearnerScenario::getSolverConfigurationFinal() {
  solver::SLESolverConfiguration solverConfigFinal;
  solverConfigFinal.eps_ = (*this)["solverFinal"]["eps"].getDouble();
  solverConfigFinal.maxIterations_ = (*this)["solverFinal"]["maxIterations"].getUInt();
  solverConfigFinal.threshold_ = (*this)["solverFinal"]["threshold"].getDouble();
  if ((*this)["solverFinal"]["type"].get().compare("CG")) {
    solverConfigFinal.type_ = solver::SLESolverType::CG;
  } else if ((*this)["solverFinal"]["type"].get().compare("BiCGSTAB")) {
    solverConfigFinal.type_ = solver::SLESolverType::BiCGSTAB;
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

/*
void LearnerScenario::writeToFile(std::string fileName) {
  std::ofstream file(fileName);

  if (file.is_open()) {
    file << "lambda=" << lambda << std::endl;
    file << "datasetFileName=" << datasetFileName << std::endl;

    file << "grid.dim="
         << "0 # inferred from dataset" << std::endl;
    file << "grid.level=" << gridConfig.level_ << std::endl;

    if (gridConfig.type_ == base::GridType::Linear) {
      file << "grid.type=Linear" << std::endl;
    } else if (gridConfig.type_ == base::GridType::ModLinear) {
      file << "grid.type=ModLinear" << std::endl;
    } else {
      throw;
    }

    file << "solverRefine.eps=" << solverConfigRefine.eps_ << std::endl;
    file << "solverRefine.maxIterations=" << solverConfigRefine.maxIterations_ << std::endl;
    file << "solverRefine.threshold=" << solverConfigRefine.threshold_ << std::endl;

    if (solverConfigRefine.type_ == solver::SLESolverType::CG) {
      file << "solverRefine.type=CG" << std::endl;
    } else if (solverConfigRefine.type_ == solver::SLESolverType::BiCGSTAB) {
      file << "solverRefine.type=BiCGSTAB" << std::endl;
    } else {
      throw;
    }

    file << "solverFinal.eps=" << solverConfigFinal.eps_ << std::endl;
    file << "solverFinal.maxIterations=" << solverConfigFinal.maxIterations_ << std::endl;
    file << "solverFinal.threshold=" << solverConfigFinal.threshold_ << std::endl;

    if (solverConfigFinal.type_ == solver::SLESolverType::CG) {
      file << "solverFinal.type=CG" << std::endl;
    } else if (solverConfigFinal.type_ == solver::SLESolverType::BiCGSTAB) {
      file << "solverFinal.type=BiCGSTAB" << std::endl;
    } else {
      throw;
    }

    file << "adapt.maxLevelType=" << adaptConfig.maxLevelType_ << std::endl;
    file << "adapt.noPoints=" << adaptConfig.noPoints_ << std::endl;
    file << "adapt.numRefinements=" << adaptConfig.numRefinements_ << std::endl;
    file << "adapt.percent=" << adaptConfig.percent_ << std::endl;
    file << "adapt.threshold=" << adaptConfig.threshold_ << std::endl;

  } else {
    throw;
  }

  file.close();
}
*/

template <class T>
T LearnerScenario::fromString(const std::string& s) {
  std::istringstream stream(s);
  T t;
  stream >> t;
  return t;
}

/*
void LearnerScenario::readFromFile(std::string fileName) {
  std::ifstream file(fileName);

  if (file.is_open()) {
    std::string line;

    uint32_t check = 0;

    while (std::getline(file, line)) {
      if (line.size() == 0) {
        continue;
      } else if (line[0] == '#') {
        continue;
      }

      std::vector<std::string> splitted;
      boost::split(splitted, line, boost::is_any_of("="));

      if (splitted.size() != 2) {
        throw;
      }
      std::string name = splitted[0];
      boost::algorithm::trim(name);

      // extract value and remove possible comments
      std::string unproccessedValue = splitted[1];
      boost::algorithm::trim(unproccessedValue);
      std::vector<std::string> splittedValue;
      boost::split(splittedValue, unproccessedValue, boost::is_any_of("#"));
      std::string value = splittedValue[0];
      boost::algorithm::trim(value);

      if (name.compare("lambda") == 0) {
        this->lambda = fromString<double>(value);
        check |= 0x00000001;
      } else if (name.compare("datasetFileName") == 0) {
        this->datasetFileName = value;
        check |= 0x00000002;
      } else if (name.compare("grid.dim") == 0) {
        if (value.compare("0") != 0) {
          throw;
        }
        gridConfig.dim_ = 0;
        check |= 0x00000004;
      } else if (name.compare("grid.level") == 0) {
        gridConfig.level_ = fromString<int>(value);
        check |= 0x00000008;
      } else if (name.compare("grid.type") == 0) {
        if (value.compare("Linear") == 0) {
          gridConfig.type_ = base::GridType::Linear;
        } else if (value.compare("ModLinear") == 0) {
          gridConfig.type_ = base::GridType::ModLinear;
        } else {
          throw;
        }
        check |= 0x00000010;
      } else if (name.compare("solverRefine.eps") == 0) {
        solverConfigRefine.eps_ = fromString<float_t>(value);
        check |= 0x00000020;
      } else if (name.compare("solverRefine.maxIterations") == 0) {
        solverConfigRefine.maxIterations_ = fromString<size_t>(value);
        check |= 0x00000040;
      } else if (name.compare("solverRefine.threshold") == 0) {
        solverConfigRefine.threshold_ = fromString<float_t>(value);
        check |= 0x00000080;
      } else if (name.compare("solverRefine.type") == 0) {
        if (value.compare("CG") == 0) {
          solverConfigRefine.type_ = solver::SLESolverType::CG;
        } else if (name.compare("BiCGSTAB") == 0) {
          solverConfigRefine.type_ = solver::SLESolverType::BiCGSTAB;
        } else {
          throw;
        }
        check |= 0x00000100;
      } else if (name.compare("solverFinal.eps") == 0) {
        solverConfigFinal.eps_ = fromString<float_t>(value);
        check |= 0x00000200;
      } else if (name.compare("solverFinal.maxIterations") == 0) {
        solverConfigFinal.maxIterations_ = fromString<size_t>(value);
        check |= 0x00000400;
      } else if (name.compare("solverFinal.threshold") == 0) {
        solverConfigFinal.threshold_ = fromString<float_t>(value);
        check |= 0x00000800;
      } else if (name.compare("solverFinal.type") == 0) {
        if (value.compare("CG") == 0) {
          solverConfigFinal.type_ = solver::SLESolverType::CG;
        } else if (name.compare("BiCGSTAB") == 0) {
          solverConfigRefine.type_ = solver::SLESolverType::BiCGSTAB;
        } else {
          throw;
        }
        check |= 0x00001000;
      } else if (name.compare("adapt.maxLevelType") == 0) {
        adaptConfig.maxLevelType_ = fromString<bool>(value);
        check |= 0x00002000;
      } else if (name.compare("adapt.noPoints") == 0) {
        adaptConfig.noPoints_ = fromString<size_t>(value);
        check |= 0x00004000;
      } else if (name.compare("adapt.numRefinements") == 0) {
        adaptConfig.numRefinements_ = fromString<size_t>(value);
        check |= 0x00008000;
      } else if (name.compare("adapt.percent") == 0) {
        adaptConfig.percent_ = fromString<float_t>(value);
        check |= 0x00010000;
      } else if (name.compare("adapt.threshold") == 0) {
        adaptConfig.threshold_ = fromString<float_t>(value);
        check |= 0x00020000;
      }
    }

    if (check != 0x0003FFFF) {
      throw;
    }

  } else {
    throw;
  }

  file.close();

  this->isInitialized = true;
}
*/

bool LearnerScenario::hasTestsetConfiguration() {
  return (*this)["testset"]["hasTestDataset"].getBool();
}

void LearnerScenario::setTestsetConfiguration(datadriven::TestsetConfiguration& testsetConfig) {
  (*this).replaceDictAttr("testset");
  (*this)["testset"].replaceIDAttr("hasTestDataset", testsetConfig.hasTestDataset);
  (*this)["testset"].replaceTextAttr("testFileName", testsetConfig.datasetFileName);
  (*this)["testset"].replaceIDAttr("expectedMSE", testsetConfig.expectedMSE);
  (*this)["testset"].replaceIDAttr("expectedLargestDifference",
                                   testsetConfig.expectedLargestDifference);
}

datadriven::TestsetConfiguration LearnerScenario::getTestsetConfiguration() {
  TestsetConfiguration testsetConfig;
  testsetConfig.hasTestDataset = (*this)["testset"]["hasTestDataset"].getBool();
  testsetConfig.datasetFileName = (*this)["testset"]["testFileName"].get();
  testsetConfig.expectedMSE = (*this)["testset"]["expectedMSE"].getDouble();
  testsetConfig.expectedLargestDifference =
      (*this)["testset"]["expectedLargestDifference"].getDouble();
  return testsetConfig;
}

}  // namespace datadriven
}  // namespace SGPP
