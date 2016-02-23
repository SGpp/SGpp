// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <string>
#include <vector>

#include "sgpp/datadriven/application/LearnerScenario.hpp"

namespace SGPP {
namespace datadriven {

// TODO(pfandedd): should be changed to json

LearnerScenario::LearnerScenario() : isInitialized(false), lambda(0.0), datasetFileName("") {}

LearnerScenario::LearnerScenario(std::string scenarioFileName) : isInitialized(true) {
  this->readFromFile(scenarioFileName);
}

LearnerScenario::LearnerScenario(std::string datasetFileName, double lambda,
                                 SGPP::base::RegularGridConfiguration gridConfig,
                                 SGPP::solver::SLESolverConfiguration SLESolverConfigRefine,
                                 SGPP::solver::SLESolverConfiguration SLESolverConfigFinal,
                                 SGPP::base::AdpativityConfiguration adaptConfig)
    : isInitialized(true),
      lambda(lambda),
      datasetFileName(datasetFileName),
      gridConfig(gridConfig),
      solverConfigRefine(SLESolverConfigRefine),
      solverConfigFinal(SLESolverConfigFinal),
      adaptConfig(adaptConfig) {}

std::string LearnerScenario::getDatasetFileName() { return this->datasetFileName; }

double LearnerScenario::getLambda() { return this->lambda; }

SGPP::base::RegularGridConfiguration LearnerScenario::getGridConfig() { return gridConfig; }

SGPP::solver::SLESolverConfiguration LearnerScenario::getSolverConfigurationRefine() {
  return solverConfigRefine;
}

SGPP::solver::SLESolverConfiguration LearnerScenario::getSolverConfigurationFinal() {
  return solverConfigFinal;
}

SGPP::base::AdpativityConfiguration LearnerScenario::getAdaptivityConfiguration() {
  return adaptConfig;
}

void LearnerScenario::writeToFile(std::string fileName) {
  std::ofstream file(fileName);

  if (file.is_open()) {
    file << "lambda=" << lambda << std::endl;
    file << "datasetFileName=" << datasetFileName << std::endl;

    file << "grid.dim="
         << "0 # inferred from dataset" << std::endl;
    file << "grid.level=" << gridConfig.level_ << std::endl;

    if (gridConfig.type_ == SGPP::base::GridType::Linear) {
      file << "grid.type=Linear" << std::endl;
    } else if (gridConfig.type_ == SGPP::base::GridType::ModLinear) {
      file << "grid.type=ModLinear" << std::endl;
    } else {
      throw;
    }

    file << "solverRefine.eps=" << solverConfigRefine.eps_ << std::endl;
    file << "solverRefine.maxIterations=" << solverConfigRefine.maxIterations_ << std::endl;
    file << "solverRefine.threshold=" << solverConfigRefine.threshold_ << std::endl;

    if (solverConfigRefine.type_ == SGPP::solver::SLESolverType::CG) {
      file << "solverRefine.type=CG" << std::endl;
    } else if (solverConfigRefine.type_ == SGPP::solver::SLESolverType::BiCGSTAB) {
      file << "solverRefine.type=BiCGSTAB" << std::endl;
    } else {
      throw;
    }

    file << "solverFinal.eps=" << solverConfigFinal.eps_ << std::endl;
    file << "solverFinal.maxIterations=" << solverConfigFinal.maxIterations_ << std::endl;
    file << "solverFinal.threshold=" << solverConfigFinal.threshold_ << std::endl;

    if (solverConfigFinal.type_ == SGPP::solver::SLESolverType::CG) {
      file << "solverFinal.type=CG" << std::endl;
    } else if (solverConfigFinal.type_ == SGPP::solver::SLESolverType::BiCGSTAB) {
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

template <class T>
T LearnerScenario::fromString(const std::string& s) {
  std::istringstream stream(s);
  T t;
  stream >> t;
  return t;
}

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
          gridConfig.type_ = SGPP::base::GridType::Linear;
        } else if (value.compare("ModLinear") == 0) {
          gridConfig.type_ = SGPP::base::GridType::ModLinear;
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
          solverConfigRefine.type_ = SGPP::solver::SLESolverType::CG;
        } else if (name.compare("BiCGSTAB") == 0) {
          solverConfigRefine.type_ = SGPP::solver::SLESolverType::BiCGSTAB;
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
          solverConfigFinal.type_ = SGPP::solver::SLESolverType::CG;
        } else if (name.compare("BiCGSTAB") == 0) {
          solverConfigRefine.type_ = SGPP::solver::SLESolverType::BiCGSTAB;
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
}  // namespace datadriven
}  // namespace SGPP
