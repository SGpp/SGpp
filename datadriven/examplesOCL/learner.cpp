// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/MetaLearner.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/globaldef.hpp>

#include <boost/program_options.hpp>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

void validate(boost::any& v, const std::vector<std::string>& values,
              sgpp::base::GridType* target_type, int) {
  // Make sure no previous assignment to 'a' was made.
  boost::program_options::validators::check_first_occurrence(v);
  // Extract the first string from 'values'. If there is more than
  // one string, it's an error, and exception will be thrown.
  const std::string& s = boost::program_options::validators::get_single_string(values);

  if (s.compare("Linear") == 0) {
    v = sgpp::base::GridType::Linear;
  } else if (s.compare("ModLinear") == 0) {
    v = sgpp::base::GridType::ModLinear;
  } else {
    throw boost::program_options::validation_error(
        boost::program_options::validation_error::invalid_option_value);
  }
}
}  // namespace base
}  // namespace sgpp

namespace sgpp {
namespace solver {

void validate(boost::any& v, const std::vector<std::string>& values,
              sgpp::solver::SLESolverType* target_type, int) {
  // Make sure no previous assignment to 'a' was made.
  boost::program_options::validators::check_first_occurrence(v);
  // Extract the first string from 'values'. If there is more than
  // one string, it's an error, and exception will be thrown.
  const std::string& s = boost::program_options::validators::get_single_string(values);

  if (s.compare("CG") == 0) {
    v = sgpp::solver::SLESolverType::CG;
  } else if (s.compare("BiCGSTAB") == 0) {
    v = sgpp::solver::SLESolverType::BiCGSTAB;
  } else {
    throw boost::program_options::validation_error(
        boost::program_options::validation_error::invalid_option_value);
  }
}
}  // namespace solver
}  // namespace sgpp

namespace sgpp {
namespace datadriven {

void validate(boost::any& v, const std::vector<std::string>& values,
              sgpp::datadriven::OperationMultipleEvalType* target_type, int) {
  // Make sure no previous assignment to 'a' was made.
  boost::program_options::validators::check_first_occurrence(v);
  // Extract the first string from 'values'. If there is more than
  // one string, it's an error, and exception will be thrown.
  const std::string& s = boost::program_options::validators::get_single_string(values);

  if (s.compare("DEFAULT") == 0) {
    v = sgpp::datadriven::OperationMultipleEvalType::DEFAULT;
  } else if (s.compare("STREAMING") == 0) {
    v = sgpp::datadriven::OperationMultipleEvalType::STREAMING;
  } else if (s.compare("SUBSPACELINEAR") == 0) {
    v = sgpp::datadriven::OperationMultipleEvalType::SUBSPACELINEAR;
  } else {
    throw boost::program_options::validation_error(
        boost::program_options::validation_error::invalid_option_value);
  }
}

void validate(boost::any& v, const std::vector<std::string>& values,
              sgpp::datadriven::OperationMultipleEvalSubType* target_type, int) {
  // Make sure no previous assignment to 'a' was made.
  boost::program_options::validators::check_first_occurrence(v);
  // Extract the first string from 'values'. If there is more than
  // one string, it's an error, and exception will be thrown.
  const std::string& s = boost::program_options::validators::get_single_string(values);

  if (s.compare("DEFAULT") == 0) {
    v = sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT;
  } else if (s.compare("COMBINED") == 0) {
    v = sgpp::datadriven::OperationMultipleEvalSubType::COMBINED;
  } else if (s.compare("OCL") == 0) {
    v = sgpp::datadriven::OperationMultipleEvalSubType::OCL;
  } else if (s.compare("OCLMASK") == 0) {
    v = sgpp::datadriven::OperationMultipleEvalSubType::OCLMASKMP;
  } else if (s.compare("OCLMP") == 0) {
    v = sgpp::datadriven::OperationMultipleEvalSubType::OCLMP;
  } else if (s.compare("SIMPLE") == 0) {
    v = sgpp::datadriven::OperationMultipleEvalSubType::SIMPLE;
  } else if (s.compare("OCLFASTMULTIPLATFORM") == 0) {
    v = sgpp::datadriven::OperationMultipleEvalSubType::OCLFASTMP;
  } else {
    throw boost::program_options::validation_error(
        boost::program_options::validation_error::invalid_option_value);
  }
}

enum class LearnerMode { LEARN, LEARNCOMPARE, LEARNTEST };

void validate(boost::any& v, const std::vector<std::string>& values,
              sgpp::datadriven::LearnerMode* target_type, int) {
  // Make sure no previous assignment to 'a' was made.
  boost::program_options::validators::check_first_occurrence(v);
  // Extract the first string from 'values'. If there is more than
  // one string, it's an error, and exception will be thrown.
  const std::string& s = boost::program_options::validators::get_single_string(values);

  if (s.compare("LEARN") == 0) {
    v = sgpp::datadriven::LearnerMode::LEARN;
  } else if (s.compare("LEARNCOMPARE") == 0) {
    v = sgpp::datadriven::LearnerMode::LEARNCOMPARE;
  } else if (s.compare("LEARNTEST") == 0) {
    v = sgpp::datadriven::LearnerMode::LEARNTEST;
  } else {
    throw boost::program_options::validation_error(
        boost::program_options::validation_error::invalid_option_value);
  }
}
}  // namespace datadriven
}  // namespace sgpp

int main(int argc, char* argv[]) {
  // std::string fileName = "debugging.arff";
  std::string trainingFileName = "DR5_train.arff";
  //  std::string fileName = "friedman2_90000.arff";
  //  std::string fileName = "bigger.arff";

  std::string testFileName = "";

  sgpp::datadriven::LearnerMode learnerMode = sgpp::datadriven::LearnerMode::LEARN;

  // only relevant for LEARNTEST-mode
  bool isRegression = true;

  double lambda = 0.000001;

  bool verbose = true;

  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::solver::SLESolverConfiguration SLESolverConfigRefine;
  sgpp::solver::SLESolverConfiguration SLESolverConfigFinal;
  sgpp::base::AdaptivityConfiguration adaptivityConfig;

  // setup grid
  gridConfig.dim_ = 0;    // dim is inferred from the data
  gridConfig.level_ = 7;  // base level
  gridConfig.type_ = sgpp::base::GridType::Linear;

  // setup adaptivity
  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.numRefinementPoints_ = 80;
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.percent_ = 200.0;
  adaptivityConfig.refinementThreshold_ = 0.0;

  // setup solver during refinement
  SLESolverConfigRefine.eps_ = 0;
  SLESolverConfigRefine.maxIterations_ = 5;
  SLESolverConfigRefine.threshold_ = -1.0;
  SLESolverConfigRefine.type_ = sgpp::solver::SLESolverType::CG;

  // setup solver for final step
  SLESolverConfigFinal.eps_ = 0;
  SLESolverConfigFinal.maxIterations_ = 5;
  SLESolverConfigFinal.threshold_ = -1.0;
  SLESolverConfigFinal.type_ = sgpp::solver::SLESolverType::CG;

  // operation type
  sgpp::datadriven::OperationMultipleEvalType type =
      sgpp::datadriven::OperationMultipleEvalType::DEFAULT;
  sgpp::datadriven::OperationMultipleEvalSubType subType =
      sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT;

  // parse command line options
  boost::program_options::options_description description("Allowed options");
  description.add_options()

      // general options
      ("help", "display help")("trainingFileName",
                               boost::program_options::value<std::string>(&trainingFileName),
                               "training data set as an arff file")(
          "testFileName", boost::program_options::value<std::string>(&testFileName),
          "test dataset as an arff file (only for LEARNTEST-mode)")(
          "isRegression", boost::program_options::value<bool>(&isRegression),
          "true for regression, false for classification, default is "
          "regression, only relevant for LEARNTEST- and LEARN-mode")(
          "lambda", boost::program_options::value<double>(&lambda),
          "regularization parameter for learning")(
          "verbose", boost::program_options::value<bool>(&verbose), "do verbose learning")(
          "learnerMode", boost::program_options::value<sgpp::datadriven::LearnerMode>(&learnerMode),
          "mode of operation: LEARN -> only learning (for performance tests), "
          "LEARNTEST -> learn and use a test dataset, LEARNCOMPARE -> learn "
          "and compare with a reference implementation and compare evaluations")

      // grid setup options
      ("grid.level", boost::program_options::value<int>(&gridConfig.level_),
       "level of the initial regular grid")(
          "grid.type", boost::program_options::value<sgpp::base::GridType>(&gridConfig.type_),
          "type of the grid to be used")

      // adaptivity options
      //    ("adaptivity.maxLevelType",
      //    boost::program_options::value<bool>(&adaptivityConfig.maxLevelType_),
      //            "DON'T KNOW WHAT THIS IS FOR")//TODO: seems to be unused,
      //            remove?
      ("adaptivityConfig.noPoints",
       boost::program_options::value<size_t>(&adaptivityConfig.numRefinementPoints_),
       "number of points to refine")(
          "adaptivityConfig.numRefinements",
          boost::program_options::value<size_t>(&adaptivityConfig.numRefinements_),
          "number of refinement steps")(
          "adaptivityConfig.percent",
          boost::program_options::value<double>(&adaptivityConfig.percent_),
          "maximum number of grid points in percent of the size of the grid "
          "that are considered for refinement")(
          "adaptivityConfig.threshold",
          boost::program_options::value<double>(&adaptivityConfig.refinementThreshold_),
          "minimum surplus value for a grid point to be considered for "
          "refinement")

      // options for the solver during refinement
      ("solverRefine.eps", boost::program_options::value<double>(&SLESolverConfigRefine.eps_),
       "error for early aborting training (set to 0 to disable)")(
          "solverRefine.maxIterations",
          boost::program_options::value<size_t>(&SLESolverConfigRefine.maxIterations_),
          "maximum number of iterations before the training is stopped")(
          "solverRefine.threshold",
          boost::program_options::value<double>(&SLESolverConfigRefine.threshold_),
          "early abort solver if this residual threshold is reached")(
          "solverRefine.type",
          boost::program_options::value<sgpp::solver::SLESolverType>(&SLESolverConfigRefine.type_),
          "the kind of solver to use")

      // options for the solver in the final step
      ("solverFinal.eps", boost::program_options::value<double>(&SLESolverConfigFinal.eps_),
       "error for early aborting training (set to 0 to disable)")(
          "solverFinal.maxIterations",
          boost::program_options::value<size_t>(&SLESolverConfigFinal.maxIterations_),
          "maximum number of iterations before the training is stopped")(
          "solverFinal.threshold",
          boost::program_options::value<double>(&SLESolverConfigRefine.threshold_),
          "early abort solver if this residual threshold is reached")(
          "solverFinal.type",
          boost::program_options::value<sgpp::solver::SLESolverType>(&SLESolverConfigFinal.type_),
          "the kind of solver to use")

      // options for the implementation type
      ("operation.type",
       boost::program_options::value<sgpp::datadriven::OperationMultipleEvalType>(&type),
       "implementation type of the operation")(
          "operation.subType",
          boost::program_options::value<sgpp::datadriven::OperationMultipleEvalSubType>(&subType),
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

  if (learnerMode == sgpp::datadriven::LearnerMode::LEARNCOMPARE) {
    std::ifstream testFile(testFileName);
    if (!testFile.good()) {
      std::cout << "test file not found" << std::endl;
      return 1;
    }
  }

  sgpp::datadriven::MetaLearner learner(gridConfig, SLESolverConfigRefine, SLESolverConfigFinal,
                                        adaptivityConfig, lambda, verbose);

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(type, subType);

  if (learnerMode == sgpp::datadriven::LearnerMode::LEARN) {
    // only execute learning (no comparisons or tests, for performance
    // measurements)
    learner.learn(configuration, trainingFileName, isRegression);
  } else if (learnerMode == sgpp::datadriven::LearnerMode::LEARNCOMPARE) {
    // execute learning with the specified configuration and use the
    // implementation from base as comparison
    // result grids are compared by sampling the domain (again with a grid) and
    // comparing the evaluated values
    learner.learnAndCompare(configuration, trainingFileName, 8);
  } else if (learnerMode == sgpp::datadriven::LearnerMode::LEARNTEST) {
    // test the learned function with a test dataset (no cross-validation yet)
    learner.learnAndTest(configuration, trainingFileName, testFileName, isRegression);
  }

  return 0;
}
