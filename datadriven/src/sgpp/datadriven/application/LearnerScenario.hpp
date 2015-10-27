/*
 * LearnerScenario.hpp
 *
 *  Created on: Oct 27, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/solver/TypesSolver.hpp>

namespace SGPP {
namespace datadriven {

class LearnerScenario {
private:
    bool isInitialized;

    //variables for the scenario
    size_t level;SGPP::base::GridType gridType;
    std::string datasetFileName;
    sg::base::RegularGridConfiguration gridConfig;
    sg::solver::SLESolverConfiguration SLESolverConfigRefine;
    sg::solver::SLESolverConfiguration SLESolverConfigFinal;
    sg::base::AdpativityConfiguration adaptConfig;
public:

    LearnerScenario() :
            isInitialized(false), level(0), gridType(SGPP::base::Linear), datasetFileName("") {
    }

    LearnerScenario(size_t level,
    SGPP::base::GridType gridType, std::string datasetFileName, sg::base::RegularGridConfiguration gridConfig,
            sg::solver::SLESolverConfiguration SLESolverConfigRefine,
            sg::solver::SLESolverConfiguration SLESolverConfigFinal, sg::base::AdpativityConfiguration adaptConfig) :
            isInitialized(true), level(level), gridType(gridType), datasetFileName(datasetFileName), gridConfig(
                    gridConfig), SLESolverConfigRefine(SLESolverConfigRefine), SLESolverConfigFinal(
                    SLESolverConfigFinal), adaptConfig(adaptConfig) {

    }

    size_t getLevel() {
        return this->level;
    }

    std::string getDatasetFileName() {
        return this->datasetFileName;
    }

    sg::base::RegularGridConfiguration getGridConfig() {
        return gridConfig;
    }

    sg::solver::SLESolverConfiguration getSLESolverConfigurationRefine() {
        return SLESolverConfigRefine;
    }

    sg::solver::SLESolverConfiguration getSLESolverConfigurationFinal() {
        return SLESolverConfigFinal;
    }

    sg::base::AdpativityConfiguration getAdaptivityConfiguration() {
        return adaptConfig;
    }

};

}
}
