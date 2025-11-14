// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// Example file demonstrating the use of IterativeGridGeneratorRitterNovak
// and IterativeGridGeneratorFullAdaptiveRitterNovak from the SGpp optimization
// module.

// To fully enable this example, ensure that the following dependencies are installed (otherwise you
// will get a slightly limited functionality):
// 1. [libgp](https://github.com/mblum/libgp#building-and-testing)
// 2. [BayesOpt](https://rmcantin.github.io/bayesopt/html/install.html)
// Additionally, compile SGpp with the `USE_LIBGP` and `USE_BAYESOPT` options enabled.
// Note: Before compiling SGpp, verify that libgp and BayesOpt are installed system-wide.

// define if a bayesopt is available for hyperparameter optimization
#ifdef USE_BAYESOPT
#include <bayesopt/bayesopt.hpp>
#endif

// Other library headers
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/sle/solver/Auto.hpp>
#include <sgpp/base/tools/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorFullAdaptiveRitterNovak.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Rastrigin.hpp>

// C system headers
#include <cassert>
#include <cmath>
#include <ctime>

// C++ system headers
#include <algorithm>  // For std::sort, std::remove, std::max
#include <array>      // For std::array
#include <chrono>     // For std::chrono
#include <fstream>    // For std::ofstream
#include <iomanip>    // For std::put_time
#include <iostream>   // For std::cout, std::endl, std::cerr
#include <limits>     // For std::numeric_limits
#include <memory>     // For std::unique_ptr, std::make_unique
#include <random>     // For std::mt19937, std::uniform_real_distribution
#include <sstream>    // For std::ostringstream
#include <stdexcept>  // For std::runtime_error, std::invalid_argument
#include <string>     // For std::string
#include <tuple>      // For std::make_tuple, std::get
#include <utility>    // For std::pair, std::make_pair
#include <vector>     // For std::vector

using sgpp::base::DataVector;
using sgpp::optimization::IGGFARNHelper::Hyperparameters;
using sgpp::optimization::IGGFARNHelper::operator<<;
using sgpp::optimization::IGGFARNHelper::UncertaintyBound;

#ifdef USE_BAYESOPT
/**
 * @brief A wrapper class to adapt a std::function to the bayesopt::ContinuousModel interface.
 *
 * This allows using bayesopt to optimize a function provided as a C++ lambda or std::function.
 */
class HyperParameterFunction : public bayesopt::ContinuousModel {
 private:
  std::function<double(const vectord&)> evaluateSampleFunc;

 public:
  /**
   * @brief Constructor for the wrapper.
   * @param num_hyperparameter Number of hyperparameters to optimize.
   * @param par bayesopt parameters.
   * @param evalFunc The std::function to be optimized.
   */
  HyperParameterFunction(const size_t num_hyperparameter, const bayesopt::Parameters& par,
                         const std::function<double(const vectord&)>& evalFunc)
      : ContinuousModel(num_hyperparameter, par), evaluateSampleFunc(evalFunc) {}

  /**
   * @brief Evaluates the wrapped function for a given sample.
   * @param xin The hyperparameter vector sample.
   * @return The result of the wrapped function.
   */
  double evaluateSample(const vectord& xin) override { return evaluateSampleFunc(xin); }
};

/**
 * @brief Overload for streaming bayesopt::Parameters to an output stream.
 * @param os The output stream.
 * @param parameters The bayesopt parameters to output.
 * @return The output stream.
 */
std::ostream& operator<<(std::ostream& os, const bayesopt::Parameters& parameters) {
  os << "n_init_samples: " << parameters.n_init_samples
     << ", n_iterations: " << parameters.n_iterations << ", random_seed: " << parameters.random_seed
     << ", noise: " << parameters.noise << ", crit_name: " << parameters.crit_name;
  return os;
}

/**
 * @brief Recursively transforms a vector of parameters (from a binary tree) into a vector of
 * weights.
 *
 * The input vector `p` of size N is assumed to represent N internal nodes of a
 * complete binary tree, defining split probabilities. The output `weights` vector
 * of size N+1 will contain the probabilities of reaching the N+1 leaf nodes.
 *
 * @param p           Input vector (e.g., from bayesopt) of size N.
 * @param weights     Output vector of size N+1 to store calculated weights.
 * @param current_value The probability weight accumulated so far down the tree.
 * @param root        The current node index in the conceptual binary tree.
 */
void transformWeights(const vectord& p, std::vector<double>& weights, const double current_value,
                      const size_t root) {
  if (root < p.size()) {
    transformWeights(p, weights, p(root) * current_value, 2 * root + 1);
    transformWeights(p, weights, (1 - p(root)) * current_value, 2 * root + 2);
  } else {
    weights[root - p.size()] = current_value;
  }
}

/**
 * @brief Executes a shell command and returns its standard output.
 * @param cmd The command string to execute.
 * @return The standard output of the command as a string.
 * @throws std::runtime_error if popen() fails.
 *
 * reference:
 * https://stackoverflow.com/questions/478898/how-do-i-execute-a-command-and-get-the-output-of-the-command-within-c-using-po
 */
std::string exec(const std::string& cmd) {
  std::array<char, 128> buffer;
  std::string result;
  // Use std::unique_ptr for automatic closing of the pipe
  std::unique_ptr<FILE, void (*)(FILE*)> pipe(popen(cmd.c_str(), "r"),
                                              [](FILE* f) -> void { std::ignore = pclose(f); });
  if (!pipe) {
    throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe.get()) != nullptr) {
    result += buffer.data();
  }
  return result;
}

/**
 * @brief Gathers project and environment information (time, git hash).
 * @return A string containing project information.
 */
std::string project_info() {
  std::ostringstream oss;

  const auto now = std::chrono::system_clock::now();
  const std::time_t time = std::chrono::system_clock::to_time_t(now);
  oss << "current time: " << std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S") << " - ";

  // Get current git branch
  std::string gitInfo = exec("git rev-parse --abbrev-ref HEAD");
  gitInfo.erase(std::remove(gitInfo.begin(), gitInfo.end(), '\n'), gitInfo.cend());
  oss << "git: " << gitInfo << " - ";
  // Get current git commit hash
  std::string gitHash = exec("git rev-parse --short HEAD");
  gitHash.erase(std::remove(gitHash.begin(), gitHash.end(), '\n'), gitHash.cend());
  oss << gitHash;

  std::string temp = oss.str();
  // Remove any extra newlines from the final string
  temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());
  return temp;
}

/**
 * @brief Calculates the probability of failure (pf) and mismatch probability (pm)
 * between surrogate and original function values.
 *
 * "Failure" is defined as the surrogate value being less than zero.
 * "Mismatch" is defined as the surrogate and original values having different signs.
 *
 * @param surrogate_values Vector of values from the surrogate model.
 * @param original_values  Vector of values from the true function.
 * @return A pair of UncertaintyBound objects: (pf, pm).
 */
std::pair<UncertaintyBound, UncertaintyBound> measure_error(
    const std::vector<double>& surrogate_values, const std::vector<double>& original_values) {
  const size_t samples = surrogate_values.size();
  assert(surrogate_values.size() == original_values.size());

  size_t pf_count = 0;
  size_t pm_count = 0;

  for (size_t i = 0; i < samples; ++i) {
    const double surrogate_value = surrogate_values[i];
    const double original_value = original_values[i];

    // Count as failure if surrogate < 0
    pf_count += static_cast<size_t>(surrogate_value < 0);

    // Count as mismatch if signs are different (XOR)
    pm_count += static_cast<size_t>((surrogate_value < 0) ^ (original_value < 0));
  }

  const UncertaintyBound pf = UncertaintyBound(pf_count, samples);
  const UncertaintyBound pm = UncertaintyBound(pm_count, samples);

  return std::make_pair(pf, pm);
}

/**
 * @brief Performs a single run of the adaptive grid generation and error measurement.
 *
 * This function sets up a sparse grid, configures the
 * IterativeGridGeneratorFullAdaptiveRitterNovak with Gaussian Process-based
 * adaptation, runs the grid generation, and finally measures the error
 * of the resulting surrogate model against a set of test points.
 *
 * @param points                  Test points for error evaluation.
 * @param real_function_values    True function values at the test points.
 * @param hyperparameters         Hyperparameters for the grid generator.
 * @param dimension               Problem dimension.
 * @param bSplineDegree           B-spline degree for the grid.
 * @param maxGridPoints           Maximum number of grid points to generate.
 * @return A pair of UncertaintyBound objects: (pf, pm) from measure_error.
 */
std::pair<UncertaintyBound, UncertaintyBound> single_run(
    const std::vector<sgpp::base::DataVector>& points,
    const std::vector<double>& real_function_values, const Hyperparameters& hyperparameters,
    const size_t dimension, const size_t bSplineDegree, const size_t maxGridPoints) {
  std::unique_ptr<sgpp::base::ScalarFunction> function =
      std::make_unique<sgpp::optimization::test_problems::RastriginObjective>(dimension);

  std::unique_ptr<sgpp::base::Grid> grid(
      sgpp::base::Grid::createBsplineGrid(dimension, bSplineDegree));

  std::vector<std::pair<double, double>> domain(dimension);
  for (size_t t = 0; t < dimension; ++t) {
    domain[t] = std::make_pair(0.0, 1.0);
  }
  sgpp::optimization::IterativeGridGeneratorFullAdaptiveRitterNovak gridGen(
      *function, *grid, maxGridPoints, hyperparameters, domain);

  if (!gridGen.generate()) {
    throw std::runtime_error("Grid generation failed, exiting");
  }

  const sgpp::base::DataVector alpha = gridGen.getAlpha();

  std::unique_ptr<sgpp::base::OperationEval> opEval(
      sgpp::op_factory::createOperationEvalNaive(*grid));
  // Create a lambda function for surrogate model evaluation
  const std::function<double(const sgpp::base::DataVector&)> surrogate_lambda =
      [&alpha, &opEval](const sgpp::base::DataVector& coordinates) {
        return opEval->eval(alpha, coordinates);
      };

  const size_t samples = points.size();
  std::vector<double> surrogate_values(samples, 0.0);
  // Evaluate the surrogate model at all test points
  for (size_t i = 0; i < samples; ++i) {
    surrogate_values[i] = surrogate_lambda(points[i]);
  }

  return measure_error(surrogate_values, real_function_values);
}

/**
 * @brief Calculates a threshold for a given function, such that exactly a
 * number of points corresponding to pf is under the threshold.
 *
 * @param points The set of points to evaluate.
 * @param function The function to evaluate.
 * @param pf The target probability of failure (fraction of points below threshold).
 * @return A pair containing the calculated threshold value and the actual UncertaintyBound (pf)
 * achieved.
 *
 * @throws std::invalid_argument if pf is not in [0, 1].
 */
std::pair<double, UncertaintyBound> inversePfcalculation(
    const std::vector<sgpp::base::DataVector>& points,
    const std::function<double(const sgpp::base::DataVector&)> function, const double pf) {
  assert(0 < points.size());
  if (pf < 0 || 1 < pf) {
    throw std::invalid_argument("function_pf must be in [0, 1]");
  }

  // Calculate the index corresponding to the pf quantile
  const size_t pf_count =
      static_cast<size_t>(std::round(pf * (static_cast<double>(points.size()) - 1)));

  std::vector<double> values;
  values.reserve(points.size());
  for (const sgpp::base::DataVector& coordinate : points) {
    values.push_back(function(coordinate));
  }
  // Sort values to find the quantile
  std::sort(values.begin(), values.end());

  // The threshold is the value at the pf_count index
  return std::make_pair(values[pf_count], UncertaintyBound(pf_count, points.size()));
}

/**
 * @brief Overload for streaming std::chrono::nanoseconds in a human-readable
 * (days, hours, minutes, seconds) format.
 * @param os    The output stream.
 * @param nanos The duration in nanoseconds.
 * @return The output stream.
 */
std::ostream& operator<<(std::ostream& os, const std::chrono::nanoseconds& nanos) {
  const size_t total_seconds = std::chrono::duration_cast<std::chrono::seconds>(nanos).count();
  const size_t seconds = total_seconds % 60;
  const size_t total_minutes = total_seconds / 60;
  const size_t minutes = total_minutes % 60;
  const size_t total_hours = total_minutes / 60;
  const size_t hours = total_hours % 24;
  const size_t days = total_hours / 24;

  os << days << "d " << hours << "h " << minutes << "m " << seconds << "s";
  return os;
}

/**
 * @brief Main function to perform hyperparameter optimization using bayesopt.
 *
 * This function sets up a test problem (Rastrigin), defines a set of
 * target failure probabilities, and then uses bayesopt to find the
 * optimal Hyperparameters for the IterativeGridGeneratorFullAdaptiveRitterNovak
 * that minimize the average relative mismatch error (prm) across all
 * test functions.
 */
void calculateOptimalHyperparameters() {
  const size_t dimension = 2;
  const size_t bSplineDegree = 3;
  const size_t maxGridPoints = 50;
  const size_t seed = 2406;
  // 5 hyperparameters will be optimized, resulting in 6 leaf weights
  const size_t number_hyperparameters = 5;
  const size_t maxSamples = 100000;

  bayesopt::Parameters parameters;
  parameters.n_init_samples = 50;
  parameters.n_iterations = 150;
  parameters.random_seed = seed;
  parameters.noise = 1e-10;
  // Use a combination of criteria for robust optimization
  parameters.crit_name = "cHedge(cEI,cLCB,cExpReturn,cOptimisticSampling)";

  const auto start = std::chrono::system_clock::now();

  std::ofstream protocol{"hyperparameter_optimisation.protocol"};
  protocol << project_info() << std::endl;

  std::unique_ptr<sgpp::base::ScalarFunction> real_function =
      std::make_unique<sgpp::optimization::test_problems::RastriginObjective>(2);

  const std::function<double(const sgpp::base::DataVector&)> real_function_lambda =
      [&](const sgpp::base::DataVector& x) { return real_function->eval(x); };

  // Generate a large set of random points for testing
  std::vector<sgpp::base::DataVector> points(maxSamples, sgpp::base::DataVector(dimension));
  std::mt19937 gen(seed);
  std::uniform_real_distribution<> dis(0.0, 1.0);
  for (size_t i = 0; i < maxSamples; ++i) {
    for (size_t d = 0; d < dimension; ++d) {
      points[i][d] = dis(gen);
    }
  }

  // Pre-calculate real function values at these points
  std::vector<double> real_function_values(maxSamples, 0.0);
  for (size_t i = 0; i < maxSamples; ++i) {
    real_function_values[i] = real_function_lambda(points[i]);
  }

  // Create a set of test functions with different failure probability thresholds
  std::vector<std::pair<double, UncertaintyBound>> test_functions{};
  for (size_t i = 1; i <= 7; ++i) {
    const double pf = 1.0 / static_cast<double>(1 << i);  // 1/2, 1/4, ..., 1/128
    test_functions.push_back(inversePfcalculation(points, real_function_lambda, pf));
  }

  std::cout << "run hyperparameter optimisation for " << test_functions.size() << " functions"
            << std::endl;

  double current_best_value = std::numeric_limits<double>::infinity();
  Hyperparameters current_best_gamma{};

  // This is the objective function for bayesopt
  const std::function<double(const vectord&)>& evalFunc = [&](const vectord& p_old) {
    // p_old has 5 values (split probabilities)
    std::vector<double> p(p_old.size() + 1);  // p will have 6 values (leaf weights)
    transformWeights(p_old, p, 1.0, 0);

    const Hyperparameters hyperparameters{p[0], p[1], p[2], p[3], p[4], p[5]};

    std::vector<std::tuple<UncertaintyBound, UncertaintyBound, UncertaintyBound>> results(
        test_functions.size());
    // Run single_run for each test function in parallel
#pragma omp parallel for
    for (size_t i = 0; i < test_functions.size(); ++i) {
      std::pair<UncertaintyBound, UncertaintyBound> pf_pm_result = single_run(
          points, real_function_values, hyperparameters, dimension, bSplineDegree, maxGridPoints);
      const auto& surrogate_pf = pf_pm_result.first;
      const auto& surrogate_pm = pf_pm_result.second;
      // Calculate relative mismatch probability
      const auto& surrogate_prm = surrogate_pm.div(test_functions[i].second);

      results[i] = std::make_tuple(surrogate_pf, surrogate_pm, surrogate_prm);
    }

    // This counter is static, it persists across calls to evalFunc
    // It is safe because it's outside the parallel loop
    static size_t counter = 0;
    protocol << "# counter: " << counter++ << " / "
             << (parameters.n_init_samples + parameters.n_iterations) << '\n';
    protocol << hyperparameters << '\n';

    // Calculate the average relative mismatch error
    double total_error = 0;
    size_t entries = 0;
    for (size_t i = 0; i < test_functions.size(); ++i) {
      const auto& result_tuple = results[i];
      const auto& pf = std::get<0>(result_tuple);
      const auto& pm = std::get<1>(result_tuple);
      const auto& prm = std::get<2>(result_tuple);

      total_error += prm.value();
      ++entries;

      protocol << test_functions[i].first << ", " << test_functions[i].second << ": surrogate: ("
               << pf << ", " << pm << ", " << prm << ")\n";
    }
    total_error /= static_cast<double>(entries);

    // Track the best hyperparameters found so far
    if (total_error < current_best_value) {
      current_best_value = total_error;
      current_best_gamma = hyperparameters;
    }
    protocol << "total_error: " << 100 * total_error
             << "% (current best: " << 100 * current_best_value << "%, " << current_best_gamma
             << ")\n"
             << std::endl;

    return total_error;
  };

  HyperParameterFunction opt(number_hyperparameters, parameters, evalFunc);

  // Set bounds for the 5 split parameters to [0, 1]
  vectord lb(number_hyperparameters);
  vectord ub(number_hyperparameters);
  for (size_t t = 0; t < number_hyperparameters; ++t) {
    lb(t) = 0.0;
    ub(t) = 1.0;
  }
  opt.setBoundingBox(lb, ub);

  opt.initializeOptimization();
  // Run the optimization loop
  for (size_t i = 0; i < parameters.n_iterations; ++i) {
    opt.stepOptimization();
  }

  const auto end = std::chrono::system_clock::now();

  protocol << project_info() << std::endl;
  protocol << "# duration: " << (end - start) << std::endl;
  protocol << "# " << parameters << std::endl;
  protocol << "best hyperparameters: " << current_best_gamma << std::endl;
  std::cout << "best hyperparameters: " << current_best_gamma << std::endl;
}
#endif

/**
 * @brief Main function to run the grid generation examples.
 * @param argc Argument count.
 * @param argv Argument values.
 * @return 0 on success, 1 on exception.
 */
int main(int argc, const char* argv[]) {
  // Silence verbose output from SG++
  sgpp::base::Printer::getInstance().setVerbosity(0);
  std::cout << "Note: For full functionality, please define USE_BAYESOPT and compile SGpp with "
               "USE_LIBGP (otherwise you can also change the code to not use gaussian processes)"
            << std::endl;
  try {
#ifdef USE_BAYESOPT
    calculateOptimalHyperparameters();
#else
    std::cout << "Note: USE_BAYESOPT not defined." << std::endl;
    std::cout << "Skipping calculateOptimalHyperparameters() example." << std::endl;
    std::cout << "This example requires the external library bayesopt." << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
#endif
  } catch (const std::exception& e) {
    std::cerr << "An exception occurred: " << e.what() << std::endl;
  }
  return 0;
}
