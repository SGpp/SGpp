#pragma once

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>

#include <set>
#include <string>

namespace sgpp {
namespace datadriven {
// Utility functions of the permutation and blow-up approach used in several classes.
namespace PermutationUtil {
/**
 * @brief Removes all elements equal to 1 from a vector and returns the obtained vector.
 *
 * @param vectorWithOnes The input vector.
 * @return std::vector<size_t>
 */
std::vector<size_t> deleteOnesFromLevelVec(std::vector<size_t> vectorWithOnes);

/**
 * @brief Permutes a vector according to a permutation suitable to transform newU into oldU.
 *
 * @param vector
 * @param oldU
 * @param newU
 * @return std::vector<size_t>
 */
std::vector<size_t> permuteVector(std::vector<size_t> vector, std::vector<size_t> currentGridLevel,
                                  std::vector<size_t> desiredGridLevel);

/**
 * @brief Returns a grid configuration with level vector without elements equal to 1 and adjusted
 * dimension.
 *
 * @param gridConfig Grid condiguration to normalize.
 * @return sgpp::base::GeneralGridConfiguration
 */
sgpp::base::GeneralGridConfiguration getNormalizedConfig(
    sgpp::base::GeneralGridConfiguration gridConfig);

/**
 * @brief Cecks whether vec1 is permutation of vec2.
 *
 * @param vec1 First input vector.
 * @param vec2 Second input vector.
 * @return true
 * @return false
 */
bool isPermutation(std::vector<size_t> vec1, std::vector<size_t> vec2);
}  // namespace PermutationUtil

class DBMatOfflinePermutable : public DBMatOffline {
 public:
  /**
   * @brief Set that contains all permutable decomposition types. If the approach is implemented for
   * a new decomposition, it has to be added to this set.
   * The set is used to check whether the permutaion and blow-up approach is applicable for a
   * configured decomposition type.
   *
   */
  static const std::set<MatrixDecompositionType> PermutableDecompositions;

  explicit DBMatOfflinePermutable(const std::string& fileName);

  /**
   * @brief Applies permutation and blow-up approach to undecomposed left-hand side matrix.
   * Meant to be used for testing primarily.
   *
   * @param baseGridConfig
   * @param desiredGridConfig
   */
  void permuteLhsMatrix(const sgpp::base::GeneralGridConfiguration& baseGridConfig,
                        const sgpp::base::GeneralGridConfiguration& desiredGridConfig);

  /**
   * @brief Applies the permutation and blow-up approach to the decomposition to match the desired
   * grid configuration.
   *
   * @param baseGridConfig Grid configuration of the current state of the offline object.
   * @param desiredGridConfig Grid configuration of the desired offline object.
   */
  virtual void permuteDecomposition(
      const sgpp::base::GeneralGridConfiguration& baseGridConfig,
      const sgpp::base::GeneralGridConfiguration& desiredGridConfig) = 0;

 protected:
  // Default constructor
  DBMatOfflinePermutable();
  /**
   * @brief Precomputation of the factors for the mapping from level and index vectors to matrix.
   * indices.
   *
   * @param level Level vector of the base grid.
   * @return std::vector<size_t>
   */
  std::vector<size_t> preComputeMatrixIndexForPoint(std::vector<size_t> level);

  /**
   * @brief Computes a permutation to transfer the desired level vector into the base level vector.
   *
   * @param baseLevel Level vector of the base object's grid.
   * @param desiredLevel Level vector of the desired object's grid.
   * @return std::vector<size_t>
   */
  std::vector<size_t> computePermutation(std::vector<size_t> baseLevel,
                                         std::vector<size_t> desiredLevel);

  /**
   * @brief Applies a permutation stored in an array to a vector.
   *
   * @param vector The vector that is to be permuted.
   * @param permutation The permutation stored in an array.
   * @return std::vector<size_t>
   */
  inline std::vector<size_t> applyPermutation(std::vector<size_t> vector,
                                              std::vector<size_t> permutation);

  /**
   * @brief Implementation of the inverse numbering function. Maps a point given by its level and
   * index vector to it's corresponding matrix index in the system matrix corresponding to the grid
   * with level vector gridLevel.
   *
   * @param level Level vector of the point
   * @param index Index vector of the point
   * @param gridLevel Level vector of the grid
   * @param preComputations Precomputed factors
   * @return size_t
   */
  size_t getMatrixIndexForPoint(std::vector<size_t> level, std::vector<size_t> index,
                                std::vector<size_t> gridLevel,
                                const std::vector<size_t>& preComputations);
  /**
   * @brief Permutes a given sytem matrix with corresponding grid config baseGridConfig to the
   * system matrix corresponding to desiredGridConfig.
   * It can be specified whether rows or columns are permuted.
   *
   * @param baseGridConfig Grid configuration of the base object
   * @param desiredGridConfig Grid configuration of the desired object
   * @param baseMatrix Base system matrix
   * @param permutedMatrix Matrix to return the system matrix of the desired object
   * @param permuteRowsOrColumns Flag to specify wheter rows or columns are permuted. For row
   * permutation set to true, for columns permutation set to false.
   */
  void permuteMatrix(const sgpp::base::GeneralGridConfiguration& baseGridConfig,
                     const sgpp::base::GeneralGridConfiguration& desiredGridConfig,
                     const sgpp::base::DataMatrix& baseMatrix,
                     sgpp::base::DataMatrix& permutedMatrix, bool permuteRowsOrColumns);

  /**
   * @brief In place application of the blow-up method to a base system matrix.
   * It can be specified wheter the base matrix is inverse or not. For inverse matrices, the inverse
   * blow-up factor is multiplied.
   *
   * @param baseGridConfig Grid configuration of the base object
   * @param desiredGridConfig   Grid configuration of the desired object
   * @param baseMatrix Base system matrix
   * @param matrixIsInverse Flag to specify whether the base matrix is inverse. Set to true if the
   * matrix is inverse. Set to false by default.
   */
  void dimensionBlowUp(const sgpp::base::GeneralGridConfiguration& baseGridConfig,
                       const sgpp::base::GeneralGridConfiguration& desiredGridConfig,
                       sgpp::base::DataMatrix& baseMatrix, bool matrixIsInverse = false);
};

}  // namespace datadriven
}  // namespace sgpp