/*
 * learner.cpp
 *
 *  Created on: Nov 2, 2015
 *      Author: pfandedd
 */

#include "boost/program_options.hpp"
//#include "boost/math/constants/constants.hpp"
//#include "boost/foreach.hpp"

#include "sgpp/datadriven/application/MetaLearner.hpp"
#include "sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp"

#include <sgpp/globaldef.hpp>

//using namespace boost::program_options;

namespace SGPP {
namespace base {

void validate(boost::any& v, const std::vector<std::string>& values,
SGPP::base::GridType* target_type, int) {
    // Make sure no previous assignment to 'a' was made.
    boost::program_options::validators::check_first_occurrence(v);
    // Extract the first string from 'values'. If there is more than
    // one string, it's an error, and exception will be thrown.
    const std::string& s = boost::program_options::validators::get_single_string(values);

    if (s.compare("Linear") == 0) {
        v = SGPP::base::GridType::Linear;
    } else if (s.compare("ModLinear") == 0) {
        v = SGPP::base::GridType::ModLinear;
    } else {
        throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
    }
}

}
}

namespace SGPP {
namespace solver {

void validate(boost::any& v, const std::vector<std::string>& values,
SGPP::solver::SLESolverType* target_type, int) {
    // Make sure no previous assignment to 'a' was made.
    boost::program_options::validators::check_first_occurrence(v);
    // Extract the first string from 'values'. If there is more than
    // one string, it's an error, and exception will be thrown.
    const std::string& s = boost::program_options::validators::get_single_string(values);

    if (s.compare("CG") == 0) {
        v = SGPP::solver::SLESolverType::CG;
    } else if (s.compare("BiCGSTAB") == 0) {
        v = SGPP::solver::SLESolverType::BiCGSTAB;
    } else {
        throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
    }
}

}
}

namespace SGPP {
namespace datadriven {

void validate(boost::any& v, const std::vector<std::string>& values,
SGPP::datadriven::OperationMultipleEvalType* target_type, int) {
    // Make sure no previous assignment to 'a' was made.
    boost::program_options::validators::check_first_occurrence(v);
    // Extract the first string from 'values'. If there is more than
    // one string, it's an error, and exception will be thrown.
    const std::string& s = boost::program_options::validators::get_single_string(values);

    if (s.compare("DEFAULT") == 0) {
        v = SGPP::datadriven::OperationMultipleEvalType::DEFAULT;
    } else if (s.compare("STREAMING") == 0) {
        v = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    } else if (s.compare("SUBSPACELINEAR") == 0) {
        v = SGPP::datadriven::OperationMultipleEvalType::SUBSPACELINEAR;
    } else {
        throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
    }
}

void validate(boost::any& v, const std::vector<std::string>& values,
SGPP::datadriven::OperationMultipleEvalSubType* target_type, int) {
    // Make sure no previous assignment to 'a' was made.
    boost::program_options::validators::check_first_occurrence(v);
    // Extract the first string from 'values'. If there is more than
    // one string, it's an error, and exception will be thrown.
    const std::string& s = boost::program_options::validators::get_single_string(values);

    if (s.compare("DEFAULT") == 0) {
        v = SGPP::datadriven::OperationMultipleEvalSubType::DEFAULT;
    } else if (s.compare("COMBINED") == 0) {
        v = SGPP::datadriven::OperationMultipleEvalSubType::COMBINED;
    } else if (s.compare("OCL") == 0) {
        v = SGPP::datadriven::OperationMultipleEvalSubType::OCL;
    } else if (s.compare("OCLMASK") == 0) {
        v = SGPP::datadriven::OperationMultipleEvalSubType::OCLMASK;
    } else if (s.compare("OCLMP") == 0) {
        v = SGPP::datadriven::OperationMultipleEvalSubType::OCLMP;
    } else if (s.compare("SIMPLE") == 0) {
        v = SGPP::datadriven::OperationMultipleEvalSubType::SIMPLE;
    } else if (s.compare("OCLFAST") == 0) {
        v = SGPP::datadriven::OperationMultipleEvalSubType::OCLFAST;
    } else if (s.compare("OCLFASTMULTIPLATFORM") == 0) {
        v = SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM;
    } else {
        throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
    }
}

enum class LearnerMode {
    LEARN, LEARNCOMPARE, LEARNTEST
};

void validate(boost::any& v, const std::vector<std::string>& values,
SGPP::datadriven::LearnerMode* target_type, int) {
    // Make sure no previous assignment to 'a' was made.
    boost::program_options::validators::check_first_occurrence(v);
    // Extract the first string from 'values'. If there is more than
    // one string, it's an error, and exception will be thrown.
    const std::string& s = boost::program_options::validators::get_single_string(values);

    if (s.compare("LEARN") == 0) {
        v = SGPP::datadriven::LearnerMode::LEARN;
    } else if (s.compare("LEARNCOMPARE") == 0) {
        v = SGPP::datadriven::LearnerMode::LEARNCOMPARE;
    } else if (s.compare("LEARNTEST") == 0) {
        v = SGPP::datadriven::LearnerMode::LEARNTEST;
    } else {
        throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
    }
}

}
}

int main(int argc, char *argv[]) {

//std::string fileName = "debugging.arff";
    std::string trainingFileName = "DR5_train.arff";
//  std::string fileName = "friedman2_90000.arff";
//  std::string fileName = "bigger.arff";

    std::string testFileName = "";

/*    SGPP::datadriven::LearnerMode learnerMode = SGPP::datadriven::LearnerMode::LEARN;

    //only relevant for LEARNTEST-mode
    bool isRegression = true;

    double lambda = 0.000001;

    bool verbose = true;

    sg::base::RegularGridConfiguration gridConfig;
    sg::solver::SLESolverConfiguration SLESolverConfigRefine;
    sg::solver::SLESolverConfiguration SLESolverConfigFinal;
    sg::base::AdpativityConfiguration adaptConfig;

// setup grid
    gridConfig.dim_ = 0;    //dim is inferred from the data
    gridConfig.level_ = 7;    //base level
    gridConfig.type_ = SGPP::base::GridType::Linear;

// setup adaptivity
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

// setup solver during refinement
    SLESolverConfigRefine.eps_ = 0;
    SLESolverConfigRefine.maxIterations_ = 5;
    SLESolverConfigRefine.threshold_ = -1.0;
    SLESolverConfigRefine.type_ = SGPP::solver::SLESolverType::CG;

// setup solver for final step
    SLESolverConfigFinal.eps_ = 0;
    SLESolverConfigFinal.maxIterations_ = 5;
    SLESolverConfigFinal.threshold_ = -1.0;
    SLESolverConfigFinal.type_ = SGPP::solver::SLESolverType::CG;

// operation type
    SGPP::datadriven::OperationMultipleEvalType type = SGPP::datadriven::OperationMultipleEvalType::DEFAULT;
    SGPP::datadriven::OperationMultipleEvalSubType subType = SGPP::datadriven::OperationMultipleEvalSubType::DEFAULT;

//parse command line options
    boost::program_options::options_description description("Allowed options");
    description.add_options()

    //general options
    ("help", "display help")("trainingFileName", boost::program_options::value<std::string>(&trainingFileName),
            "training data set as an arff file")("testFileName",
            boost::program_options::value<std::string>(&testFileName),
            "test dataset as an arff file (only for LEARNTEST-mode)")
            ("isRegression", boost::program_options::value<bool>(&isRegression),
            "true for regression, false for classification, default is regression, only relevant for LEARNTEST- and LEARN-mode")("lambda",
            boost::program_options::value<SGPP::float_t>(&lambda), "regularization parameter for learning")("verbose",
            boost::program_options::value<bool>(&verbose), "do verbose learning")("learnerMode",
            boost::program_options::value<SGPP::datadriven::LearnerMode>(&learnerMode),
            "mode of operation: LEARN -> only learning (for performance tests), LEARNTEST -> learn and use a test dataset, LEARNCOMPARE -> learn and compare with a reference implementation and compare evaluations")

    //grid setup options
    ("grid.level", boost::program_options::value<int>(&gridConfig.level_), "level of the initial regular grid")(
            "grid.type", boost::program_options::value<SGPP::base::GridType>(&gridConfig.type_),
            "type of the grid to be used")

    //adaptivity options
//    ("adaptivity.maxLevelType", boost::program_options::value<bool>(&adaptConfig.maxLevelType_),
//            "DON'T KNOW WHAT THIS IS FOR")//TODO: seems to be unused, remove?
    ("adaptConfig.noPoints", boost::program_options::value<size_t>(&adaptConfig.noPoints_),
            "number of points to refine")("adaptConfig.numRefinements",
            boost::program_options::value<size_t>(&adaptConfig.numRefinements_), "number of refinement steps")(
            "adaptConfig.percent", boost::program_options::value<SGPP::float_t>(&adaptConfig.percent_),
            "maximum number of grid points in percent of the size of the grid that are considered for refinement")(
            "adaptConfig.threshold", boost::program_options::value<SGPP::float_t>(&adaptConfig.threshold_),
            "minimum surplus value for a grid point to be considered for refinement")

    //options for the solver during refinement
    ("solverRefine.eps", boost::program_options::value<SGPP::float_t>(&SLESolverConfigRefine.eps_),
            "error for early aborting training (set to 0 to disable)")("solverRefine.maxIterations",
            boost::program_options::value<size_t>(&SLESolverConfigRefine.maxIterations_),
            "maximum number of iterations before the training is stopped")("solverRefine.threshold",
            boost::program_options::value<SGPP::float_t>(&SLESolverConfigRefine.threshold_),
            "early abort solver if this residual threshold is reached")
    ("solverRefine.type", boost::program_options::value<SGPP::solver::SLESolverType>(&SLESolverConfigRefine.type_),
            "the kind of solver to use")

    //options for the solver in the final step
    ("solverFinal.eps", boost::program_options::value<SGPP::float_t>(&SLESolverConfigFinal.eps_),
            "error for early aborting training (set to 0 to disable)")("solverFinal.maxIterations",
            boost::program_options::value<size_t>(&SLESolverConfigFinal.maxIterations_),
            "maximum number of iterations before the training is stopped")("solverFinal.threshold",
            boost::program_options::value<SGPP::float_t>(&SLESolverConfigRefine.threshold_),
            "early abort solver if this residual threshold is reached")
    ("solverFinal.type", boost::program_options::value<SGPP::solver::SLESolverType>(&SLESolverConfigFinal.type_),
            "the kind of solver to use")

    //options for the implementation type
    ("operation.type", boost::program_options::value<SGPP::datadriven::OperationMultipleEvalType>(&type),
            "implementation type of the operation")("operation.subType",
            boost::program_options::value<SGPP::datadriven::OperationMultipleEvalSubType>(&subType),
            "implementation sub type of the operation");

    boost::program_options::variables_map variables_map;

    boost::program_options::parsed_options options = parse_command_line(argc, argv, description);
    boost::program_options::store(options, variables_map);
    boost::program_options::notify(variables_map);

    if (variables_map.count("help")) {
      std::cout << description << std::endl;
      return 0;
    }

    std::ifstream trainingFile(trainingFileName);
    if (!trainingFile.good()) {
        std::cout << "training file not found" << std::endl;
        return 1;
    }

    if (learnerMode == SGPP::datadriven::LearnerMode::LEARNCOMPARE) {
        std::ifstream testFile(testFileName);
        if (!testFile.good()) {
            std::cout << "test file not found" << std::endl;
            return 1;
        }
    }

    SGPP::datadriven::MetaLearner learner(gridConfig, SLESolverConfigRefine, SLESolverConfigFinal, adaptConfig, lambda,
            verbose);

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(type, subType);

    if (learnerMode == SGPP::datadriven::LearnerMode::LEARN) {
// only execute learning (no comparisons or tests, for performance measurements)
        learner.learn(configuration, trainingFileName, isRegression);
    } else if (learnerMode == SGPP::datadriven::LearnerMode::LEARNCOMPARE) {
// execute learning with the specified configuration and use the implementation from base as comparison
// result grids are compared by sampling the domain (again with a grid) and comparing the evaluated values
        learner.learnAndCompare(configuration, trainingFileName, 8);
    } else if (learnerMode == SGPP::datadriven::LearnerMode::LEARNTEST) {
// test the learned function with a test dataset (no cross-validation yet)
        learner.learnAndTest(configuration, trainingFileName, testFileName, isRegression);
		}*/

    return 0;
}

