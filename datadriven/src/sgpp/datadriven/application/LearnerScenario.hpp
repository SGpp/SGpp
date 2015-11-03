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
//    size_t level;SGPP::base::GridType gridType;
    double lambda;
    std::string datasetFileName;
    sg::base::RegularGridConfiguration gridConfig;
    sg::solver::SLESolverConfiguration solverConfigRefine;
    sg::solver::SLESolverConfiguration solverConfigFinal;
    sg::base::AdpativityConfiguration adaptConfig;
public:

    LearnerScenario() :
            isInitialized(false),
            lambda(0.0),
//        level(0), gridType(SGPP::base::Linear),
            datasetFileName("") {
    }

    LearnerScenario(
//            size_t level,
//    SGPP::base::GridType gridType,
    std::string datasetFileName, double lambda, SGPP::base::RegularGridConfiguration gridConfig,
    SGPP::solver::SLESolverConfiguration SLESolverConfigRefine,
    SGPP::solver::SLESolverConfiguration SLESolverConfigFinal, SGPP::base::AdpativityConfiguration adaptConfig) :
            isInitialized(true), lambda(lambda),
//            level(level), gridType(gridType),
            datasetFileName(datasetFileName), gridConfig(
                    gridConfig), solverConfigRefine(SLESolverConfigRefine), solverConfigFinal(
                    SLESolverConfigFinal), adaptConfig(adaptConfig) {

    }

//    size_t getLevel() {
//        return this->gridConfig.level_;
//    }

    std::string getDatasetFileName() {
        return this->datasetFileName;
    }

    double getLambda() {
        return this->lambda;
    }

    sg::base::RegularGridConfiguration getGridConfig() {
        return gridConfig;
    }

    sg::solver::SLESolverConfiguration getSolverConfigurationRefine() {
        return solverConfigRefine;
    }

    sg::solver::SLESolverConfiguration getSolverConfigurationFinal() {
        return solverConfigFinal;
    }

    sg::base::AdpativityConfiguration getAdaptivityConfiguration() {
        return adaptConfig;
    }

};

}
}
