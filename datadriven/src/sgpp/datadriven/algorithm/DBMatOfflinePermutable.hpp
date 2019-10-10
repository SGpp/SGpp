#pragma once

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>

#include <set>
#include <string>

namespace sgpp {
namespace datadriven {
namespace PermutationUtil {
std::vector<size_t> deleteOnesFromLevelVec(std::vector<size_t> vectorWithOnes);
std::vector<size_t> permutateVector(std::vector<size_t> vector, std::vector<size_t> oldU,
                                    std::vector<size_t> newU);
sgpp::base::GeneralGridConfiguration getNormalizedConfig(
    sgpp::base::GeneralGridConfiguration gridConfig);

bool isPermutation(std::vector<size_t> vec1, std::vector<size_t> vec2);
}  // namespace PermutationUtil

class DBMatOfflinePermutable : public DBMatOffline {
 public:
  /**
   * Constant that contains all decompositio, for which the permutation method is currently
   * implemented.
   */
  static const std::set<MatrixDecompositionType> PermutableDecompositions;

  explicit DBMatOfflinePermutable(const std::string& fileName);
  /**
   * Permutates the base matrix and multiplies the dimension blowup factor to obtain the desired
   * system matrix.
   */
  void permutateLhsMatrix(sgpp::base::GeneralGridConfiguration baseGridConfig,
                          sgpp::base::GeneralGridConfiguration desiredGridCOnfig);

  virtual void permutateDecomposition(sgpp::base::GeneralGridConfiguration baseGridConfig,
                                      sgpp::base::GeneralGridConfiguration desiredGridCOnfig) = 0;

 protected:
  DBMatOfflinePermutable();

  std::vector<size_t> preComputeMatrixIndexForPoint(std::vector<size_t> level);

  std::vector<size_t> computePermutation(std::vector<size_t> baseLevel,
                                                            std::vector<size_t> desiredLevel);
  inline std::vector<size_t> applyPermutation(std::vector<size_t> vector,
                                              std::vector<size_t> permutation);

  size_t getMatrixIndexForPoint(std::vector<size_t> level, std::vector<size_t> index,
                                std::vector<size_t> gridLevel,
                                const std::vector<size_t>& preComputations);

  void permutateMatrix(sgpp::base::GeneralGridConfiguration baseGridConfig,
                       sgpp::base::GeneralGridConfiguration desiredGridCOnfig,
                       const sgpp::base::DataMatrix& baseMatrix,
                       sgpp::base::DataMatrix& permutatedMatrix, bool permutateRowsOrColumns);

  void dimensionBlowUp(sgpp::base::GeneralGridConfiguration baseGridConfig,
                       sgpp::base::GeneralGridConfiguration desiredGridCOnfig,
                       sgpp::base::DataMatrix& baseMatrix, bool matrixIsInverse = false);
};

}  // namespace datadriven
}  // namespace sgpp