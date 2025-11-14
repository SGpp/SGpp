// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/gridgen/IterativeGridGeneratorFullAdaptiveRitterNovak.hpp>

#include <sgpp/base/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/base/function/scalar/InterpolantScalarFunctionGradient.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/sle/solver/IterativeGaussianElimination.hpp>
#include <sgpp/base/tools/sle/system/HierarchisationSLE.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/gridgen/HashRefinementMultiple.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iomanip>  // For std::setprecision
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef USE_LIBGP
#include <gp.h>
#include <rprop.h>
#endif

namespace sgpp {
namespace optimization {

namespace {
// Type aliases for clarity
using LeftRightPoint = std::array<base::GridPoint, 2>;
using CriterionRankSingle = std::vector<std::vector<double>>;
using CriterionRankLR = std::vector<std::vector<std::array<double, 2>>>;
using CriterionRanks = std::pair<CriterionRankSingle, CriterionRankLR>;

// Confidence level for Wilson score interval
constexpr double accuracy = 0.999;

/**
 * @brief CDF of the standard normal distribution.
 * @param z The z-score.
 * @return P(Z <= z).
 */
double phi(const double z) { return 0.5 * (1.0 + std::erf(z / std::sqrt(2.0))); }

/**
 * @brief Computes the quantile of the standard normal distribution (inverse of the CDF).
 * @param p Probability (must be in (0, 1)).
 * @return The corresponding z-score.
 */
static double inverse_phi(const double p) {
  if (p <= 0.0) {
    return -std::numeric_limits<double>::infinity();
  } else if (p >= 1.0) {
    return std::numeric_limits<double>::infinity();
  }

  const double tolerance_p = 1e-12;

  double low_z = -10.0;           // = inverse_phi(tolerance_p)
  double high_z = 10.0;           // = inverse_phi(1 - tolerance_p)
  const int max_iterations = 50;  // static_cast<int>(std::ceil(std::log2((high_z - low_z) /
                                  // tolerance_p))) + 5; // 5 safety puffer

  double mid_z;
  for (int i = 0; i < max_iterations; ++i) {
    mid_z = low_z + (high_z - low_z) / 2.0;  // prevent overflow

    const double diff_p = phi(mid_z) - p;
    if (std::abs(diff_p) < tolerance_p) {
      return mid_z;
    }

    if (diff_p < 0) {
      low_z = mid_z;
    } else {
      high_z = mid_z;
    }
  }

  throw std::runtime_error(
      "inverse_phi: Failed to converge to a solution within the maximum number of iterations.");
}

/**
 * @brief Calculates the Wilson score confidence interval for a binomial proportion.
 * It is more reliable than the Wald interval, especially for p near 0 or 1.
 *
 * @param p     Observed proportion.
 * @param N     Total number of samples.
 * @param gamma Confidence level (e.g., 0.999).
 * @return A pair containing the lower and upper bounds of the confidence interval.
 */
static std::pair<double, double> confidence_interval_wilson(double p, size_t N, double gamma) {
  assert(0.0 <= p && p <= 1.0);
  assert(0.0 <= gamma && gamma <= 1.0);
  assert(0 < N);

  const double alpha = 1.0 - gamma;
  const double c = inverse_phi(1.0 - alpha / 2.0);
  const double cSquared = c * c;

  const double denominator = 2.0 * (static_cast<double>(N) + cSquared);
  const double numerator_center = 2.0 * static_cast<double>(N) * p + cSquared;
  const double numerator_margin =
      c * std::sqrt(cSquared + 4.0 * static_cast<double>(N) * p * (1.0 - p));

  const double lower = (numerator_center - numerator_margin) / denominator;
  const double upper = (numerator_center + numerator_margin) / denominator;

  return std::make_pair(std::max(0.0, lower), std::min(1.0, upper));
}

/**
 * @brief Creates an empty rank matrix for single-value criteria.
 * @param dimension The grid dimension.
 * @param currentN The current number of grid points.
 * @return A 2D vector of the specified size.
 */
CriterionRankSingle emptyCriterionRankSingle(const size_t dimension, const size_t currentN) {
  return std::vector<std::vector<double>>(dimension, std::vector<double>(currentN));
}

/**
 * @brief Creates an empty rank matrix for left/right criteria.
 * @param dimension The grid dimension.
 * @param currentN The current number of grid points.
 * @return A 2D vector of arrays of the specified size.
 */
CriterionRankLR emptyCriterionRankLR(const size_t dimension, const size_t currentN) {
  return std::vector<std::vector<std::array<double, 2>>>(
      dimension, std::vector<std::array<double, 2>>(currentN));
}

/**
 * @brief Calculates the 1-based rank of each element in a vector of values.
 * Rank is defined as #{j | values[j] <= values[i]}.
 * Handles ties by assigning sequential ranks based on sorting order.
 *
 * @param values A vector of double values.
 * @return A vector of ranks corresponding to the input values.
 */
std::vector<double> calculateRank(const std::vector<double>& values) {
  std::vector<size_t> p(values.size());
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(), [&](size_t i, size_t j) { return values[i] < values[j]; });

  std::vector<double> rank(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    rank[p[i]] = static_cast<double>(i + 1);  // +1 for 1-based rank
  }
  return rank;
}

/**
 * @brief Calculates a single rank per (point, dimension) pair by taking the minimum of the
 * left and right child criteria values before ranking.
 *
 * @param values Vector of [left_value, right_value] pairs.
 * @return A vector of single ranks.
 */
std::vector<double> calculateRankSingle(const std::vector<std::array<double, 2>>& values) {
  std::vector<double> stripped(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    stripped[i] = std::min(values[i][0], values[i][1]);
  }
  return calculateRank(stripped);
}

/**
 * @brief Calculates separate ranks for left and right children and scales them.
 *
 * @param values Vector of [left_value, right_value] pairs.
 * @return A vector of [left_rank, right_rank] pairs.
 */
std::vector<std::array<double, 2>> calculateRankLR(
    const std::vector<std::array<double, 2>>& values) {
  std::vector<double> stripped;
  stripped.reserve(values.size() * 2);
  for (const auto& value : values) {
    stripped.push_back(value[0]);
    stripped.push_back(value[1]);
  }
  std::vector<double> rank = calculateRank(stripped);
  std::vector<std::array<double, 2>> result(values.size());
  const double N = static_cast<double>(values.size());
  const double a = (N - 1.0) / (2.0 * N - 1.0);
  const double b = N / (2.0 * N - 1.0);
  for (size_t i = 0; i < values.size(); ++i) {
    result[i][0] = a * static_cast<double>(rank[2 * i + 0]) + b;
    result[i][1] = a * static_cast<double>(rank[2 * i + 1]) + b;
  }
  return result;
}

/**
 * @brief Converts a matrix of criterion values into rank matrices.
 * @param fX A matrix of criterion values [dimension][point][left/right].
 * @return A pair of matrices: one for single ranks, one for left/right ranks.
 */
CriterionRanks toRank(const std::vector<std::vector<std::array<double, 2>>>& fX) {
  if (fX.empty() || fX[0].empty()) {
    return std::make_pair(CriterionRankSingle(), CriterionRankLR());
  }
  CriterionRanks rank = std::make_pair(emptyCriterionRankSingle(fX.size(), fX[0].size()),
                                       emptyCriterionRankLR(fX.size(), fX[0].size()));
  for (size_t t = 0; t < fX.size(); ++t) {
    rank.first[t] = calculateRankSingle(fX[t]);
    rank.second[t] = calculateRankLR(fX[t]);
  }
  return rank;
}

/**
 * @brief Determines the grid points of the two children that would be created by refining
 * grid point `i` in dimension `t`. This is a "what-if" analysis and does not
 * permanently modify the grid.
 *
 * @param gridStorage The grid's storage.
 * @param i Index of the grid point to consider for refinement.
 * @param t Dimension to consider for refinement.
 * @return An array containing the two child grid points.
 */
LeftRightPoint getRefinableChildren(base::GridStorage& gridStorage, const size_t i,
                                    const size_t t) {
  const bool isLeaf = gridStorage.getPoint(i).isLeaf();
  gridStorage.getPoint(i).setLeaf(false);

  const size_t oldGridSize = gridStorage.getSize();

  HashRefinementMultiple refinement;
  refinement.refineGridpoint1D(gridStorage, gridStorage.getPoint(i), t);
  assert(oldGridSize + 2 == gridStorage.getSize());

  std::array<base::GridPoint, 2> addedPoints;
  addedPoints[0] = gridStorage.getPoint(oldGridSize);
  addedPoints[1] = gridStorage.getPoint(oldGridSize + 1);

  // Backtrack: remove the newly added points to not alter the grid
  while (oldGridSize < gridStorage.getSize()) {
    gridStorage.deleteLast();
  }

  gridStorage.getPoint(i).setLeaf(isLeaf);

  return addedPoints;
}

bool isSamePoint(const base::HashGridPoint& firstPoint, const base::HashGridPoint& secondPoint) {
  if (firstPoint.getDimension() != secondPoint.getDimension()) {
    return false;
  }

  for (size_t t = 0; t < firstPoint.getDimension(); ++t) {
    base::HashGridPoint::level_type firstLevel, secondLevel;
    base::HashGridPoint::index_type firstIndex, secondIndex;

    firstPoint.get(t, firstLevel, firstIndex);
    secondPoint.get(t, secondLevel, secondIndex);

    if (firstLevel != secondLevel || firstIndex != secondIndex) {
      return false;
    }
  }

  return true;
}

/**
 * @brief Updates the list of potential children for all refineable grid points.
 *
 * @param gridStorage The grid's storage.
 * @param dimension The grid dimension.
 * @param currentN The current number of grid points.
 * @param refineableGridPoints The matrix of child points to update.
 * @param ignore A matrix to flag children that should be ignored (e.g., max level exceeded).
 * @param maxLevel The maximum allowed grid level.
 * @param domain The domain of the function.
 */
void updateRefineableGridPoints(base::GridStorage& gridStorage, const size_t dimension,
                                const size_t currentN,
                                std::vector<std::vector<LeftRightPoint>>& refineableGridPoints,
                                std::vector<std::vector<std::array<bool, 2>>>& ignore,
                                base::level_t maxLevel,
                                const std::vector<std::pair<double, double>>& domain) {
  for (size_t t = 0; t < dimension; ++t) {
    const size_t oldN = refineableGridPoints[t].size();
    refineableGridPoints[t].resize(currentN);

    for (size_t i = 0; i < currentN; i++) {
      bool needUpdate = false;
      if (oldN <= i) {
        needUpdate = true;
      } else {
        for (size_t q = oldN; q < currentN && !needUpdate; ++q) {
          for (size_t lr = 0; lr < refineableGridPoints[t][i].size() && !needUpdate; ++lr) {
            needUpdate = isSamePoint(refineableGridPoints[t][i][lr], gridStorage[q]);
          }
        }
      }

      const LeftRightPoint newPoints = getRefinableChildren(gridStorage, i, t);
      for (size_t lr = 0; lr < 2; ++lr) {
        base::DataVector coordinates(dimension);
        newPoints[lr].getStandardCoordinates(coordinates);

        const bool oldIgnore = (i < oldN) ? ignore[t][i][lr] : false;
        const base::GridPoint oldPoint =
            (i < oldN) ? refineableGridPoints[t][i][lr] : base::GridPoint(dimension);

        ignore[t][i][lr] = maxLevel < newPoints[lr].getLevelMax();
        for (size_t q = 0; q < dimension; ++q) {
          if (coordinates[q] < domain[q].first || domain[q].second < coordinates[q]) {
            ignore[t][i][lr] = true;
            break;
          }
        }

        if (!ignore[t][i][lr]) {
          refineableGridPoints[t][i][lr] = newPoints[lr];
        }

        if (!needUpdate) {
          assert(oldIgnore == ignore[t][i][lr] &&
                 isSamePoint(oldPoint, refineableGridPoints[t][i][lr]));
        }
      }
    }
  }
}

/**
 * @brief Evaluates the interpolant's gradient at parent points and value/gradient at child
 * points.
 *
 * @param grid The current grid.
 * @param alpha The current surplus vector.
 * @param currentN The current number of grid points.
 * @param dimension The grid dimension.
 * @param refineableGridPoints The matrix of child points to evaluate.
 * @param ignore A matrix flagging children that should be ignore.
 * @return A pair containing (gradients at parents, {value, gradient} at children).
 */
std::pair<std::vector<base::DataVector>,
          std::vector<std::vector<std::array<std::pair<double, base::DataVector>, 2>>>>
evalValues(base::Grid& grid, const base::DataVector& alpha, const size_t currentN,
           const size_t dimension,
           const std::vector<std::vector<LeftRightPoint>>& refineableGridPoints,
           const std::vector<std::vector<std::array<bool, 2>>>& ignore) {
  const base::GridStorage& gridStorage = grid.getStorage();

  std::vector<base::DataVector> functionGradients(currentN, base::DataVector(dimension));
  std::vector<std::vector<std::array<std::pair<double, base::DataVector>, 2>>> childrenValues(
      dimension,
      std::vector<std::array<std::pair<double, base::DataVector>, 2>>(
          currentN, {{{0.0, base::DataVector(dimension)}, {0.0, base::DataVector(dimension)}}}));

  base::InterpolantScalarFunctionGradient ftGradient(grid, alpha);
  base::DataVector coordinates(dimension);
  base::DataVector gradient(dimension);

  for (size_t i = 0; i < currentN; ++i) {
    gridStorage.getPoint(i).getStandardCoordinates(coordinates);
    ftGradient.eval(coordinates, functionGradients[i]);

    for (size_t t = 0; t < dimension; ++t) {
      for (size_t lr = 0; lr < 2; ++lr) {
        if (ignore[t][i][lr]) {
          continue;
        }

        refineableGridPoints[t][i][lr].getStandardCoordinates(coordinates);
        childrenValues[t][i][lr].first = ftGradient.eval(coordinates, gradient);
        childrenValues[t][i][lr].second = gradient;
      }
    }
  }

  return std::make_pair(std::move(functionGradients), std::move(childrenValues));
}

/**
 * @brief Calculates ranks based on the sum of level vectors (L1 norm).
 * This criterion penalizes points that are already highly refined.
 *
 * @param gridStorage The grid's storage.
 * @param dimension The grid dimension.
 * @param currentN The current number of grid points.
 * @param degree The number of times each point has been refined.
 * @return A pair of rank matrices.
 */
CriterionRanks evalRankLevelDegree(const base::GridStorage& gridStorage, const size_t dimension,
                                   const size_t currentN,
                                   const std::vector<std::vector<std::array<size_t, 2>>>& degree) {
  // levelSum[i] is the 1-norm of the level vector of the i-th grid point.
  // Intuition: Don't refine points which are already highly refined (iff level is high),
  // therefore `levelSum` is independent from the refined dimension
  const bool useL2Norm = false;
  const bool useL1Norm = true;
  const bool useOneDimension = false;

  std::vector<std::vector<double>> levelSum(dimension, std::vector<double>(currentN, 0.0));
  for (size_t t = 0; t < dimension; ++t) {
    for (size_t i = 0; i < currentN; ++i) {
      const base::GridPoint& gp = gridStorage[i];

      if (useL2Norm) {
        for (size_t d = 0; d < dimension; ++d) {
          levelSum[t][i] +=
              static_cast<double>(gp.getLevel(d)) * static_cast<double>(gp.getLevel(d));
        }
        levelSum[t][i] = std::sqrt(levelSum[t][i]);
      } else if (useL1Norm) {
        for (size_t d = 0; d < dimension; ++d) {
          levelSum[t][i] += static_cast<double>(gp.getLevel(d));
        }
      } else if (useOneDimension) {
        levelSum[t][i] = static_cast<double>(gp.getLevel(t));
      } else {
        assert(false);
      }
    }
  }

  CriterionRankSingle rankSingle = emptyCriterionRankSingle(dimension, currentN);
  CriterionRankLR rankLR = emptyCriterionRankLR(dimension, currentN);
  for (size_t t = 0; t < dimension; ++t) {
    for (size_t i = 0; i < currentN; ++i) {
      for (size_t lr = 0; lr < 2; ++lr) {
        rankLR[t][i][lr] = levelSum[t][i] + static_cast<double>(degree[t][i][lr]);
      }
      rankSingle[t][i] = std::min(rankLR[t][i][0], rankLR[t][i][1]);
    }
  }
  return std::make_pair(rankSingle, rankLR);
}

/**
 * @brief Calculates ranks based on the absolute function value at the parent point.
 * This is the classic surplus-based refinement criterion.
 *
 * @param functionValues The surplus vector.
 * @param dimension The grid dimension.
 * @param currentN The current number of grid points.
 * @return A pair of rank matrices.
 */
CriterionRanks evalRankThreshold(const base::DataVector& functionValues, const size_t dimension,
                                 const size_t currentN) {
  CriterionRankLR fXThreshold = emptyCriterionRankLR(dimension, currentN);
  for (size_t t = 0; t < dimension; ++t) {
    for (size_t i = 0; i < currentN; ++i) {
      fXThreshold[t][i][0] = fXThreshold[t][i][1] = std::abs(functionValues[i]);
    }
  }
  return toRank(fXThreshold);
}

/**
 * @brief Calculates ranks based on the absolute interpolated function value at the child points.
 *
 * @param dimension The grid dimension.
 * @param currentN The current number of grid points.
 * @param childrenValues The evaluated values at child points.
 * @return A pair of rank matrices.
 */
CriterionRanks evalRankExtendedThreshold(
    const size_t dimension, size_t currentN,
    const std::vector<std::vector<std::array<std::pair<double, base::DataVector>, 2>>>&
        childrenValues) {
  CriterionRankLR fXExtendedThreshold = emptyCriterionRankLR(dimension, currentN);
  for (size_t t = 0; t < dimension; ++t) {
    for (size_t i = 0; i < currentN; ++i) {
      fXExtendedThreshold[t][i][0] = std::abs(childrenValues[t][i][0].first);
      fXExtendedThreshold[t][i][1] = std::abs(childrenValues[t][i][1].first);
    }
  }
  return toRank(fXExtendedThreshold);
}

/**
 * @brief Calculates ranks based on the gradient at the parent point.
 * The value is `1 / (|grad_t| + 1)` to prioritize refinement in directions of high gradient.
 *
 * @param dimension The grid dimension.
 * @param currentN The current number of grid points.
 * @param functionGradients The gradients at parent points.
 * @return A pair of rank matrices.
 */
CriterionRanks evalRankGradient(const size_t dimension, const size_t currentN,
                                const std::vector<base::DataVector>& functionGradients) {
  CriterionRankLR fXGradient = emptyCriterionRankLR(dimension, currentN);
  for (size_t t = 0; t < dimension; ++t) {
    for (size_t i = 0; i < currentN; i++) {
      const double value = 1.0 / (std::abs(functionGradients[i][t]) + 1.0);
      fXGradient[t][i][0] = fXGradient[t][i][1] = value;
    }
  }
  return toRank(fXGradient);
}

/**
 * @brief Calculates ranks based on the change in gradient from parent to child.
 * The value is `1 / (|grad_parent_t - grad_child_t| + 1)`.
 *
 * @param dimension The grid dimension.
 * @param currentN The current number of grid points.
 * @param functionGradients The gradients at parent points.
 * @param childrenValues The values and gradients at child points.
 * @return A pair of rank matrices.
 */
CriterionRanks evalRankGradientDeviation(
    const size_t dimension, const size_t currentN,
    const std::vector<base::DataVector>& functionGradients,
    const std::vector<std::vector<std::array<std::pair<double, base::DataVector>, 2>>>&
        childrenValues) {
  CriterionRankLR fXGradientDeviation = emptyCriterionRankLR(dimension, currentN);
  for (size_t t = 0; t < dimension; ++t) {
    for (size_t i = 0; i < currentN; ++i) {
      const double parentGradientComp = functionGradients[i][t];
      fXGradientDeviation[t][i][0] =
          1.0 / (std::abs(parentGradientComp - childrenValues[t][i][0].second[t]) + 1.0);
      fXGradientDeviation[t][i][1] =
          1.0 / (std::abs(parentGradientComp - childrenValues[t][i][1].second[t]) + 1.0);
    }
  }
  return toRank(fXGradientDeviation);
}

#ifdef USE_LIBGP
void initializeGaussianProcess(std::unique_ptr<libgp::GaussianProcess>& gp,
                               const size_t dimension) {
  // Define the covariance function (kernel)
  // libgp uses a string format:
  // - CovSum: Adds subsequent kernels.
  // - CovMatern3iso: Matern kernel with smoothness 3/2 (isotropic).
  // - CovNoise: White noise kernel.
  std::string covf_str = "CovSum( CovMatern3iso, CovNoise )";
  gp = std::make_unique<libgp::GaussianProcess>(dimension, covf_str);

  // How quickly the function changes
  const double length_scale = 0.05;
  // Overall variance of the function
  const double signal_variance = 1.0;
  // Assumed noise in the observations
  const double noise_variance = 1e-3;

  // libgp requires hyperparameters to be set in log-space.
  Eigen::VectorXd params(gp->covf().get_param_dim());
  params << std::log(length_scale), std::log(signal_variance), std::log(noise_variance);

  // These are initial guesses - ideally, they should be optimized!
  gp->covf().set_loghyper(params);

  base::Printer::getInstance().printStatusUpdate("Using Kernel: " + gp->covf().to_string());
  std::stringstream ss;
  ss << "Initial hyperparameters: [" << length_scale << ", " << signal_variance << ", "
     << noise_variance << "]";
  base::Printer::getInstance().printStatusUpdate(ss.str());
}

std::vector<std::vector<std::array<std::pair<double, double>, 2>>> evaluateGaussianProcess(
    std::unique_ptr<libgp::GaussianProcess>& gp,
    const std::vector<std::vector<std::array<base::GridPoint, 2>>>& refineableGridPoints,
    const std::vector<std::vector<std::array<bool, 2>>>& ignore, const size_t dimension) {
  // Optimize hyperparameters (e.g., using RProp)
  libgp::RProp rProp;
  rProp.init();
  // Maximize log-likelihood, 5 iterations, no verbosity
  rProp.maximize(&*gp, 5, false);

  // Assume, that the noise parameter is always at index 2, which corresponds to the order in
  // initializeGaussianProcess: (length_scale, signal_variance, noise_variance).
  // If the kernel string is changed, this index WILL be wrong.
  static const size_t kNoiseParamIndex = 2;
  const double kResetNoiseVariance = 1e-3;

  Eigen::VectorXd params{gp->covf().get_loghyper()};
  // Check if params vector is large enough, to avoid out-of-bounds access
  if (params.size() > static_cast<Eigen::Index>(kNoiseParamIndex)) {
    params[kNoiseParamIndex] = std::log(kResetNoiseVariance);
    gp->covf().set_loghyper(params);
  } else {
    throw std::runtime_error(
        "Error: Could not reset noise variance. Hyperparameter vector size is " +
        std::to_string(params.size()) + ", but index " + std::to_string(kNoiseParamIndex) +
        " was expected.");
  }

  std::vector<std::vector<std::array<std::pair<double, double>, 2>>> result(
      refineableGridPoints.size());
  base::DataVector new_point(dimension);

  for (size_t t = 0; t < refineableGridPoints.size(); ++t) {
    result[t].resize(refineableGridPoints[t].size());
    for (size_t i = 0; i < refineableGridPoints[t].size(); ++i) {
      result[t][i] = std::array<std::pair<double, double>, 2>();
      for (size_t lr = 0; lr < 2; ++lr) {
        if (ignore[t][i][lr]) {
          continue;
        }

        refineableGridPoints[t][i][lr].getStandardCoordinates(new_point);

        const double predicted_mean = gp->f(new_point.data());

        // Warning: `gp->var` can return small negative values due to numerical instabilities. Use
        // std::abs to prevent sqrt(negative) which results in NaN.
        const double predicted_variance = std::abs(gp->var(new_point.data()));
        const double predicted_stddev = std::sqrt(predicted_variance);

        result[t][i][lr] = std::make_pair(predicted_mean, predicted_stddev);
      }
    }
  }

  return result;
}

/**
 * @brief Calculates ranks based on Gaussian Process uncertainty.
 *
 * The value is computed as `phi(|mu|/sigma)`. This prioritizes points with high uncertainty or a
 * mean close to zero.
 *
 * @param dimension             The grid dimension.
 * @param currentN              The current number of grid points.
 * @param uncertaintyGridPoints The GP mean and standard deviation at child points.
 * @return A pair of rank matrices.
 */
CriterionRanks evalRankGaussianProcess(
    const size_t dimension, const size_t currentN,
    const std::vector<std::vector<std::array<std::pair<double, double>, 2>>>&
        uncertaintyGridPoints) {
  CriterionRankLR fXGpZScore = emptyCriterionRankLR(dimension, currentN);

  for (size_t t = 0; t < dimension; ++t) {
    for (size_t i = 0; i < currentN; ++i) {
      for (size_t lr = 0; lr < 2; ++lr) {
        const double mean = uncertaintyGridPoints[t][i][lr].first;
        const double stddev = uncertaintyGridPoints[t][i][lr].second;
        // If stddev is effectively zero, uncertainty is minimal. Rank it high (value 1.0) to avoid
        // selection
        fXGpZScore[t][i][lr] = (stddev < 1e-9) ? 1.0 : phi(std::abs(mean) / stddev);
      }
    }
  }

  return toRank(fXGpZScore);
}
#endif
}  // namespace

IterativeGridGeneratorFullAdaptiveRitterNovak::IterativeGridGeneratorFullAdaptiveRitterNovak(
    base::ScalarFunction& f, base::Grid& grid, const size_t N,
    const IGGFARNHelper::Hyperparameters& hyperparameters,
    const std::vector<std::pair<double, double>>& domain, const bool refineLeftOrRight,
    const base::level_t initialLevel, const base::level_t maxLevel)
    : IterativeGridGenerator(f, grid, N),
      hyperparameters(hyperparameters),
      domain(domain),
      initialLevel(initialLevel),
      maxLevel(maxLevel),
      stoppingCriterionValidationPoints(),
      refineLeftOrRight(refineLeftOrRight) {}

IterativeGridGeneratorFullAdaptiveRitterNovak::~IterativeGridGeneratorFullAdaptiveRitterNovak() {}

const IGGFARNHelper::Hyperparameters&
IterativeGridGeneratorFullAdaptiveRitterNovak::getHyperparameters() const {
  return this->hyperparameters;
}

void IterativeGridGeneratorFullAdaptiveRitterNovak::setHyperparameters(
    const IGGFARNHelper::Hyperparameters& hyperparameters) {
  this->hyperparameters = hyperparameters;
}

base::level_t IterativeGridGeneratorFullAdaptiveRitterNovak::getInitialLevel() const {
  return this->initialLevel;
}

void IterativeGridGeneratorFullAdaptiveRitterNovak::setInitialLevel(base::level_t initialLevel) {
  this->initialLevel = initialLevel;
}

base::level_t IterativeGridGeneratorFullAdaptiveRitterNovak::getMaxLevel() const {
  return this->maxLevel;
}

void IterativeGridGeneratorFullAdaptiveRitterNovak::setMaxLevel(base::level_t maxLevel) {
  this->maxLevel = maxLevel;
}

void IterativeGridGeneratorFullAdaptiveRitterNovak::activateStoppingCriterion(
    const std::vector<base::DataVector>& validationPoints) {
  if (validationPoints.empty()) {
    throw std::runtime_error("Validation points cannot be empty.");
  }
  this->stoppingCriterionValidationPoints = validationPoints;
}

void IterativeGridGeneratorFullAdaptiveRitterNovak::disableStoppingCriterion() {
  this->stoppingCriterionValidationPoints.clear();
}

bool IterativeGridGeneratorFullAdaptiveRitterNovak::isStoppingCriterionEnabled() const {
  return !this->stoppingCriterionValidationPoints.empty();
}

const base::DataVector& IterativeGridGeneratorFullAdaptiveRitterNovak::getAlpha() const {
  return this->alpha;
}

bool IterativeGridGeneratorFullAdaptiveRitterNovak::generate() {
  base::Printer::getInstance().printStatusBegin(
      "Adaptive grid generation (Full Adaptive Ritter-Novak)...");

  base::GridStorage& gridStorage = grid.getStorage();
  const size_t dim = f.getNumberOfParameters();

  // Generate initial grid
  grid.getGenerator().regular(initialLevel);
  size_t currentN = gridStorage.getSize();

  // Resize functionValues to max grid size
  const size_t maxGridSize = std::max(N, currentN);
  functionValues.resize(maxGridSize);
  functionValues.setAll(0.0);

  // degree[t][i] counts refinements of point i in dimension t
  std::vector<std::vector<std::array<size_t, 2>>> degree(
      dim, std::vector<std::array<size_t, 2>>(maxGridSize, {{0, 0}}));
  // ignore[t][i] flags points that cannot be refined further
  std::vector<std::vector<std::array<bool, 2>>> ignore(
      dim, std::vector<std::array<bool, 2>>(maxGridSize, {{false, false}}));

  // Evaluate f on initial grid points
  evalFunction();

  base::sle_solver::IterativeGaussianElimination fastMatrixSolver{};
  std::unique_ptr<IGGFARNHelper::StoppingCriterion> stoppingCriterion;
  if (this->isStoppingCriterionEnabled()) {
    stoppingCriterion =
        std::make_unique<IGGFARNHelper::StoppingCriterion>(this->stoppingCriterionValidationPoints);
  }

  size_t k = 0;  // Iteration counter
  base::HierarchisationSLE hierSLE(grid);
  std::vector<std::vector<LeftRightPoint>> refineableGridPoints(dim);

#ifdef USE_LIBGP
  std::unique_ptr<libgp::GaussianProcess> gp;
  initializeGaussianProcess(gp, dim);
  for (size_t i = 0; i < currentN; ++i) {
    gp->add_pattern(gridStorage.getPointCoordinates(i).data(), functionValues[i]);
  }
#endif

  const auto start = std::chrono::high_resolution_clock::now();
  const double epsilon = 1e-9;  // For float comparison

  while (currentN < N) {
    // Status printing
    {
      const auto current = std::chrono::high_resolution_clock::now();
      const auto elapsed =
          std::chrono::duration_cast<std::chrono::milliseconds>(current - start).count();

      std::stringstream sstream;
      sstream << std::fixed << std::setprecision(2);
      sstream << (static_cast<double>(currentN) / static_cast<double>(N) * 100.0)
              << "%\t " + std::to_string(currentN) << " / " << N
              << ",\t elapsed time: " << static_cast<double>(elapsed) / 1000.0 << "s";
      base::Printer::getInstance().printStatusUpdate(sstream.str());
    }

    // Hierarchize
    base::DataVector functionValuesCurrent(functionValues);
    functionValuesCurrent.resize(currentN);
    if (!fastMatrixSolver.iterativeSolve(hierSLE, functionValuesCurrent, this->alpha)) {
      throw std::runtime_error("Hierarchisation failed during adaptive grid generation.");
    }

    if (stoppingCriterion && stoppingCriterion->shouldStop(this->grid, this->alpha)) {
      break;
    }

    updateRefineableGridPoints(gridStorage, dim, currentN, refineableGridPoints, ignore, maxLevel,
                               this->domain);

    const auto& evalPair =
        evalValues(grid, this->alpha, currentN, dim, refineableGridPoints, ignore);
    const auto& functionGradients = evalPair.first;
    const auto& childrenValues = evalPair.second;

    // Calculate all rank criteria
    std::vector<std::pair<double, CriterionRanks>> weightedRanks;

    if (std::abs(this->hyperparameters.getLevelDegree()) > epsilon) {
      weightedRanks.push_back({this->hyperparameters.getLevelDegree(),
                               evalRankLevelDegree(gridStorage, dim, currentN, degree)});
    }
    if (std::abs(this->hyperparameters.getThreshold()) > epsilon) {
      weightedRanks.push_back({this->hyperparameters.getThreshold(),
                               evalRankThreshold(this->functionValues, dim, currentN)});
    }
    if (std::abs(this->hyperparameters.getExtendedThreshold()) > epsilon) {
      weightedRanks.push_back({this->hyperparameters.getExtendedThreshold(),
                               evalRankExtendedThreshold(dim, currentN, childrenValues)});
    }
    if (std::abs(this->hyperparameters.getGradient()) > epsilon) {
      weightedRanks.push_back({this->hyperparameters.getGradient(),
                               evalRankGradient(dim, currentN, functionGradients)});
    }
    if (std::abs(this->hyperparameters.getSecondGradient()) > epsilon) {
      weightedRanks.push_back(
          {this->hyperparameters.getSecondGradient(),
           evalRankGradientDeviation(dim, currentN, functionGradients, childrenValues)});
    }
    if (std::abs(this->hyperparameters.getGaussianProcess()) > epsilon) {
#ifdef USE_LIBGP
      const std::vector<std::vector<std::array<std::pair<double, double>, 2>>>
          uncertaintyGridPoints = evaluateGaussianProcess(gp, refineableGridPoints, ignore, dim);
      weightedRanks.push_back({this->hyperparameters.getGaussianProcess(),
                               evalRankGaussianProcess(dim, currentN, uncertaintyGridPoints)});
#else
      throw std::runtime_error(
          "Gaussian Process functionality is not enabled. Please define USE_LIBGP.");
#endif
    }

    // Determine the best refinement
    struct BestRefinement {
      size_t leftRight = 0;
      size_t index = 0;
      size_t dimension = 0;
      double value = std::numeric_limits<double>::infinity();
    } bestRefinement;

    for (size_t i = 0; i < currentN; i++) {
      for (size_t t = 0; t < dim; ++t) {
        if (this->refineLeftOrRight) {
          for (size_t lr = 0; lr < 2; ++lr) {
            if (ignore[t][i][lr]) continue;

            // Combined refinement criterion
            double g = 1.0;
            for (const auto& weightRankPair : weightedRanks) {
              g *= std::pow(static_cast<double>(weightRankPair.second.second[t][i][lr]) + 1.0,
                            weightRankPair.first);
            }

            if (g < bestRefinement.value) {
              bestRefinement.leftRight = lr, bestRefinement.index = i, bestRefinement.dimension = t,
              bestRefinement.value = g;
            }
          }
        } else {
          if (ignore[t][i][0] || ignore[t][i][1]) continue;

          // Combined refinement criterion
          double g = 1.0;
          for (const auto& weightRankPair : weightedRanks) {
            g *= std::pow(static_cast<double>(weightRankPair.second.first[t][i]) + 1.0,
                          weightRankPair.first);
          }

          if (g < bestRefinement.value) {
            bestRefinement.leftRight = 0, bestRefinement.index = i, bestRefinement.dimension = t,
            bestRefinement.value = g;
          }
        }
      }
    }

    // Refine the best point
    for (size_t lr = 0; lr < 2; ++lr) {
      if (!ignore[bestRefinement.dimension][bestRefinement.index][lr] &&
          (!this->refineLeftOrRight || bestRefinement.leftRight == lr)) {
        ++degree[bestRefinement.dimension][bestRefinement.index][lr];
        gridStorage.insert(
            refineableGridPoints[bestRefinement.dimension][bestRefinement.index][lr]);
      }
    }

    // new grid size
    const size_t newN = gridStorage.getSize();

    if (newN == currentN) {
      // size unchanged ==> point not refined (should not happen with valid settings)
      base::Printer::getInstance().printStatusEnd(
          "error: grid size unchanged, refinement failed. Stopping.");
      return false;
    } else if (newN > N) {
      // too many new points ==> undo refinement and exit
      undoRefinement(currentN);
      break;
    }

    // Evaluate function at new points
    evalFunction(currentN);

#ifdef USE_LIBGP
    for (size_t i = currentN; i < newN; ++i) {
      gp->add_pattern(gridStorage.getPointCoordinates(i).data(), functionValues[i]);
    }
#endif

    // next round
    currentN = newN;
    k++;
  }

  // delete superfluous entries in functionValues
  functionValues.resize(currentN);
  fastMatrixSolver.iterativeSolve(hierSLE, functionValues, this->alpha);

  // Final status and statistics
  std::stringstream finalMsg;
  finalMsg << "100.0% (N = " << currentN << ", k = " << k << ")\n";
  finalMsg << "Refinement counts per dimension: (";
  size_t totalRefinements = 0;
  std::vector<size_t> refinementCounter(dim, 0);
  for (size_t t = 0; t < dim; ++t) {
    for (size_t i = 0; i < currentN; ++i) {
      refinementCounter[t] += degree[t][i][0] + degree[t][i][1];
    }
    totalRefinements += refinementCounter[t];
  }
  for (size_t t = 0; t < dim; ++t) {
    finalMsg << (t == 0 ? "" : ", ") << std::fixed << std::setprecision(1)
             << (totalRefinements > 0 ? 100.0 * static_cast<double>(refinementCounter[t]) /
                                            static_cast<double>(totalRefinements)
                                      : 0.0)
             << "%";
  }
  finalMsg << ")";
  base::Printer::getInstance().printStatusEnd(finalMsg.str());

  return true;
}

namespace IGGFARNHelper {

// Definition of function from header
UncertaintyBound::UncertaintyBound(double lower, double value, double upper, double coV)
    : _lower(lower), _value(value), _upper(upper), _coV(coV) {}

UncertaintyBound::UncertaintyBound() : UncertaintyBound(0.0, 0.0, 0.0, 0.0) {}

UncertaintyBound::UncertaintyBound(const size_t samples, const size_t N)
    : _value(static_cast<double>(samples) / static_cast<double>(N)) {
  assert(samples <= N);
  assert(0 < N);

  // For a Bernoulli distribution, the coefficient of variation is sqrt((1-p)/(p*N)).
  // If p is 0, CoV is infinite. We return 1.0 as a placeholder.
  this->_coV =
      (0.0 == this->_value) ? 1.0 : std::sqrt((1.0 - _value) / (_value * static_cast<double>(N)));
  std::tie(this->_lower, this->_upper) = confidence_interval_wilson(this->_value, N, accuracy);
}

UncertaintyBound UncertaintyBound::div(const UncertaintyBound& other) const {
  const double epsilon = 1e-9;
  // Division by zero is handled by returning a safe but potentially wide interval.
  const double lower = (std::abs(other._upper) < epsilon) ? 0.0 : this->_lower / other._upper;
  const double value = (std::abs(other._value) < epsilon) ? 1.0 : this->_value / other._value;
  const double upper = (std::abs(other._lower) < epsilon) ? 1.0 : this->_upper / other._lower;

  return UncertaintyBound(lower, value, upper, -1.0);  // CoV not meaningful after division
}

double UncertaintyBound::lower() const { return _lower; }
double UncertaintyBound::value() const { return _value; }
double UncertaintyBound::upper() const { return _upper; }
double UncertaintyBound::coV() const { return _coV; }

std::ostream& operator<<(std::ostream& os, const UncertaintyBound& bounds) {
  os << "[" << 100.0 * bounds.lower() << "%, " << 100.0 * bounds.value() << "%, "
     << 100.0 * bounds.upper() << "%]";
  return os;
}

Hyperparameters::Hyperparameters(const double levelDegree, const double threshold,
                                 const double extendedThreshold, const double gradient,
                                 const double secondGradient, const double gaussianProcess)
    : levelDegree(levelDegree),
      threshold(threshold),
      extendedThreshold(extendedThreshold),
      gradient(gradient),
      secondGradient(secondGradient),
      gaussianProcess(gaussianProcess) {
  assert(0.0 <= this->levelDegree && this->levelDegree <= 1.0);
  assert(0.0 <= this->threshold && this->threshold <= 1.0);
  assert(0.0 <= this->extendedThreshold && this->extendedThreshold <= 1.0);
  assert(0.0 <= this->gradient && this->gradient <= 1.0);
  assert(0.0 <= this->secondGradient && this->secondGradient <= 1.0);
  assert(0.0 <= this->gaussianProcess && this->gaussianProcess <= 1.0);

  const double sum = this->levelDegree + this->threshold + this->extendedThreshold +
                     this->gradient + this->secondGradient + this->gaussianProcess;

  // The scaling of the hyperparameters doesn't influence the selection criterion of the next
  // point, if we scale in each step. The reason is the, that the modified RitterNovakCriterion
  // uses the hyperparameters in the exponent, such that a scaling corresponds to a constant
  // factor
  if (1e-9 < sum) {
    this->levelDegree /= sum;
    this->threshold /= sum;
    this->extendedThreshold /= sum;
    this->gradient /= sum;
    this->secondGradient /= sum;
    this->gaussianProcess /= sum;
  }
}

Hyperparameters::Hyperparameters() : Hyperparameters(0, 0, 0, 0, 0, 0) {}

double Hyperparameters::getLevelDegree() const { return this->levelDegree; }
double Hyperparameters::getThreshold() const { return this->threshold; }
double Hyperparameters::getExtendedThreshold() const { return this->extendedThreshold; }
double Hyperparameters::getGradient() const { return this->gradient; }
double Hyperparameters::getSecondGradient() const { return this->secondGradient; }
double Hyperparameters::getGaussianProcess() const { return this->gaussianProcess; }

std::ostream& operator<<(std::ostream& os, const Hyperparameters& hyperparameters) {
  std::stringstream sstream;
  sstream << std::fixed;
  sstream.precision(8);
  sstream << "(levelDegree, threshold, extendedThreshold, gradient, secondGradient, "
             "gaussianProcess) = ("
          << hyperparameters.getLevelDegree() << ", " << hyperparameters.getThreshold() << ", "
          << hyperparameters.getExtendedThreshold() << ", " << hyperparameters.getGradient() << ", "
          << hyperparameters.getSecondGradient() << ", " << hyperparameters.getGaussianProcess()
          << ")";
  os << sstream.str();
  return os;
}

StoppingCriterion::StoppingCriterion(const std::vector<base::DataVector>& monteCarloPoints)
    : monteCarloPoints(monteCarloPoints) {}

StoppingCriterion::~StoppingCriterion() {}

bool StoppingCriterion::shouldStop(base::Grid& grid, const base::DataVector& alpha) {
  // Clone the current grid and alpha for comparison
  last_grid_values.emplace_back(grid.getSize(), std::vector<double>(),
                                std::unique_ptr<base::Grid>(grid.clone()), alpha);

  // Find an older grid to compare against that has a sufficiently different number of points
  int old_grid_index = static_cast<int>(last_grid_values.size()) - 2;
  while (0 <= old_grid_index && std::get<0>(last_grid_values[static_cast<size_t>(old_grid_index)]) +
                                        stopping_criterion_count_difference >
                                    grid.getSize()) {
    --old_grid_index;
  }

  if (0 <= old_grid_index) {
    auto& new_grid_tuple = last_grid_values.back();
    base::InterpolantScalarFunction ft(*std::get<2>(new_grid_tuple), std::get<3>(new_grid_tuple));

    // The logic requires that ALL comparisons to recent-enough older grids show convergence.
    // The `&& should_stop` in the loop condition ensures this: if any comparison fails to
    // signal convergence, the loop terminates and the function returns false.
    bool should_stop = true;
    for (size_t i = static_cast<size_t>(old_grid_index);
         i < last_grid_values.size() - 1 && should_stop; ++i) {
      auto& old_grid_tuple = last_grid_values[i];
      base::InterpolantScalarFunction ft_old(*std::get<2>(old_grid_tuple),
                                             std::get<3>(old_grid_tuple));

      should_stop = false;  // Assume we don't stop, prove otherwise in this iteration
      measure_error(
          monteCarloPoints, std::get<1>(new_grid_tuple), std::get<1>(old_grid_tuple),
          [&ft](const base::DataVector& coords) { return ft.eval(coords); },
          [&ft_old](const base::DataVector& coords) { return ft_old.eval(coords); },
          [&should_stop](const UncertaintyBound& pf, const UncertaintyBound& pm,
                         size_t /*numPoints*/) -> bool {
            const UncertaintyBound approx_prm = pm.div(pf);
            if (pf.value() > 1e-9 && pm.value() > 1e-9 && pf.value() < 1.0 - 1e-9 &&
                pm.value() < 1.0 - 1e-9) {
              if (stopping_criterion_tolerance < approx_prm.lower()) {
                // Mismatch is significantly higher than tolerance -> don't stop
                return true;
              } else if (approx_prm.upper() < stopping_criterion_tolerance) {
                // Mismatch is significantly lower than tolerance -> we can stop
                should_stop = true;
                return true;
              }
            }
            return false;
          });
    }

    if (should_stop) {
      base::Printer::getInstance().printStatusUpdate("stopping criterion reached at N = " +
                                                     std::to_string(grid.getSize()));
      return true;
    }
  }
  return false;
}

std::tuple<size_t, UncertaintyBound, UncertaintyBound> StoppingCriterion::measure_error(
    const std::vector<base::DataVector>& points, std::vector<double>& surrogate_values,
    std::vector<double>& original_values,
    const std::function<double(const base::DataVector&)>& surrogate_f,
    const std::function<double(const base::DataVector&)>& original_f,
    const std::function<bool(const UncertaintyBound&, const UncertaintyBound&, size_t)>&
        early_stop_criterion) {
  const size_t batchSize = 5000;
  const size_t maxSamples = points.size();
  if (0 == maxSamples) {
    return std::make_tuple(0, UncertaintyBound(), UncertaintyBound());
  }

  UncertaintyBound pf, pm;
  size_t pf_count = 0;
  size_t pm_count = 0;
  size_t samples = 0;

  while (samples < maxSamples) {
    const size_t end_batch = std::min(samples + batchSize, maxSamples);

    // Cache values if not already computed
    const size_t surrogate_values_size_old = surrogate_values.size();
    if (surrogate_values_size_old < end_batch) {
      surrogate_values.resize(end_batch);
    }

    const size_t original_values_size_old = original_values.size();
    if (original_values_size_old < end_batch) {
      original_values.resize(end_batch);
    }

    for (size_t i = samples; i < end_batch; ++i) {
      if (surrogate_values_size_old <= i) {
        surrogate_values[i] = surrogate_f(points[i]);
      }
      if (original_values_size_old <= i) {
        original_values[i] = original_f(points[i]);
      }

      pf_count += static_cast<size_t>(surrogate_values[i] < 0.0);
      pm_count += static_cast<size_t>((surrogate_values[i] < 0.0) != (original_values[i] < 0.0));
    }

    samples = end_batch;
    pf = UncertaintyBound(pf_count, samples);
    pm = UncertaintyBound(pm_count, samples);

    if (early_stop_criterion(pf, pm, samples)) {
      break;
    }
  }
  return std::make_tuple(samples, pf, pm);
}

}  // namespace IGGFARNHelper
}  // namespace optimization
}  // namespace sgpp
