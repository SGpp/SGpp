#pragma once

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

class DBMatOfflinePermutable : public DBMatOffline {
 public:
  explicit DBMatOfflinePermutable(const std::string& fileName);
  /**
   * Permutates the base matrix and multiplies the dimension blowup factor to obtain the desired
   * system matrix.
   */
  void permutateLhsMatrix(sgpp::base::CombiGridConfiguration baseGridConfig,
                          sgpp::base::CombiGridConfiguration desiredGridCOnfig);

  virtual void permutateDecomposition(sgpp::base::CombiGridConfiguration baseGridConfig,
                                      sgpp::base::CombiGridConfiguration desiredGridCOnfig) = 0;

  DBMatOfflinePermutable();

  size_t getMatrixIndexForPoint(std::vector<size_t> level, std::vector<size_t> index,
                                std::vector<size_t> gridLevel);
  std::vector<size_t> deleteOnesFromLevelVec(std::vector<size_t> vectorWithOnes);
  std::vector<size_t> permutateVector(std::vector<size_t> vector, std::vector<size_t> oldU,
                                      std::vector<size_t> newU);

  void permutateMatrix(sgpp::base::CombiGridConfiguration baseGridConfig,
                       sgpp::base::CombiGridConfiguration desiredGridCOnfig,
                       const sgpp::base::DataMatrix& baseMatrix,
                       sgpp::base::DataMatrix& permutatedMatrix, bool permutateRowsOrColumns);

  void dimensionBlowUp(sgpp::base::CombiGridConfiguration baseGridConfig,
                       sgpp::base::CombiGridConfiguration desiredGridCOnfig,
                       sgpp::base::DataMatrix& baseMatrix, bool matrixIsInverse = false);
};

}  // namespace datadriven
}  // namespace sgpp