

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflinePermutable.hpp>

namespace sgpp {
namespace datadriven {

// Implementatio of the utility functions

std::vector<size_t> PermutationUtil::deleteOnesFromLevelVec(std::vector<size_t> vectorWithOnes) {
  std::vector<size_t> output;
  for (size_t i = 0; i < vectorWithOnes.size(); i++) {
    if (vectorWithOnes[i] != 1) output.push_back(vectorWithOnes[i]);
  }
  return output;
}

sgpp::base::GeneralGridConfiguration PermutationUtil::getNormalizedConfig(
    sgpp::base::GeneralGridConfiguration gridConfig) {
  sgpp::base::GeneralGridConfiguration normalizedConfig = gridConfig;
  std::vector<size_t> levelVecWithoutOnes =
      PermutationUtil::deleteOnesFromLevelVec(gridConfig.levelVector_);
  normalizedConfig.levelVector_.resize(levelVecWithoutOnes.size());
  normalizedConfig.levelVector_.swap(levelVecWithoutOnes);
  normalizedConfig.dim_ = levelVecWithoutOnes.size();
  return normalizedConfig;
}

std::vector<size_t> PermutationUtil::permutateVector(std::vector<size_t> vector,
                                                     std::vector<size_t> oldU,
                                                     std::vector<size_t> newU) {
  std::vector<size_t> output(vector.size(), 1);
  for (size_t i = 0; i < newU.size(); i++) {
    size_t oldIndex = -1;
    for (size_t j = 0; j < oldU.size(); j++) {
      if (oldU[j] == newU[i]) {
        oldIndex = j;
        oldU[j] = 0;
        break;
      }
    }
    if (oldIndex == -1)
      throw sgpp::base::algorithm_exception("No permuation.");
    else {
      output[i] = vector[oldIndex];
    }
  }
  return output;
}

bool PermutationUtil::isPermutation(std::vector<size_t> vec1, std::vector<size_t> vec2) {
  // only vectors with same dimesion can be permutations
  if (vec1.size() != vec2.size()) return false;

  std::vector<int> vec2_(vec1.size(), 0);
  for (size_t i = 0; i < vec2.size(); i++) {
    vec2_[i] = vec2[i];
  }

  for (auto l : vec1) {
    // iterate through base level vector to find suitable element
    for (size_t i = 0; i < vec2.size(); i++) {
      auto l_ = vec2[i];
      // if suitable elemet is found, remove it from the base vector by setting it to -1.
      if (l == l_) {
        vec2[i] = -1;
        break;
      }
      // If no suitable element has been found yet, base level vec is no permutation
      else if (i == vec2.size() - 1) {
        return false;
        break;
      }
    }
  }
  return true;
}

const std::set<MatrixDecompositionType> DBMatOfflinePermutable::PermutableDecompositions{
    MatrixDecompositionType::OrthoAdapt};

DBMatOfflinePermutable::DBMatOfflinePermutable() : DBMatOffline() {}

DBMatOfflinePermutable::DBMatOfflinePermutable(const std::string& fileName)
    : DBMatOffline(fileName) {}

std::vector<size_t> DBMatOfflinePermutable::preComputeMatrixIndexForPoint(
    std::vector<size_t> level) {
  std::vector<size_t> result(level.size() - 1);
  for (size_t i = 1; i < level.size(); i++) {
    size_t prod = 1;
    for (size_t j = 0; j < i; j++) {
      prod *= (1 << level[j]) - 1;
    }
    result[i - 1] = prod;
  }
  return result;
}

std::vector<size_t> DBMatOfflinePermutable::computePermutation(std::vector<size_t> vec1,
                                                               std::vector<size_t> vec2) {
  // result vector
  std::vector<size_t> output(vec1.size(), 0);
  for (size_t i = 0; i < vec1.size(); i++) {
    size_t oldIndex = -1;
    for (size_t j = 0; j < vec2.size(); j++) {
      if (vec2[j] == vec1[i]) {
        oldIndex = j;
        vec2[j] = 0;
        break;
      }
    }
    if (oldIndex == -1)
      throw sgpp::base::algorithm_exception("No permuation.");
    else {
      output.at(i) = oldIndex;
    }
  }
  return output;
}

inline std::vector<size_t> DBMatOfflinePermutable::applyPermutation(
    std::vector<size_t> vector, std::vector<size_t> permutation) {
  // std::assert(vector.size() == permutation.size());
  std::vector<size_t> result(vector.size());

  for (size_t i = 0; i < vector.size(); i++) {
    result[i] = vector.at(permutation.at(i));
  }
  return result;
}

size_t DBMatOfflinePermutable::getMatrixIndexForPoint(std::vector<size_t> level,
                                                      std::vector<size_t> index,
                                                      std::vector<size_t> gridLevel,
                                                      const std::vector<size_t>& preComputations) {
  if (!(level.size() == index.size() && index.size() == gridLevel.size()))
    throw sgpp::base::algorithm_exception("Vector dimensions do not match");
  // Iterate over vectors. Start with most right element unequal 1
  int lStar = -1;
  for (size_t i = 0; i < level.size(); i++) {
    if (level[i] != 1) lStar = i;
  }

  size_t result = 1;
  size_t mult = 1;
  while (lStar >= 0) {
    if (lStar == 0) {
      result += ((1 << (level[0] - 1)) - 2 + ((index[0] + 1) >> 1)) * mult;

      break;
    }
    size_t prod = preComputations.at(lStar - 1);

    result += (prod + (1 << (level[lStar] - 1)) + ((index[lStar] + 1) >> 1) - 3) * mult;

    mult *= (1 << gridLevel[lStar]) - 2;

    // determine new lStar
    int lStar_ = -1;
    for (size_t i = lStar - 1; i >= 0; i--) {
      if (level[i] > 1) {
        lStar_ = i;
        break;
      }
    }
    lStar = lStar_;
  }
  return result;
}

void DBMatOfflinePermutable::permuteMatrix(
    const sgpp::base::GeneralGridConfiguration& baseGridConfig,
    const sgpp::base::GeneralGridConfiguration& desiredGridCOnfig,
    const sgpp::base::DataMatrix& baseMatrix, sgpp::base::DataMatrix& permutatedMatrix,
    bool permutateRowsOrColums) {
  if (baseMatrix.getNrows() != permutatedMatrix.getNrows() ||
      baseMatrix.getNcols() != permutatedMatrix.getNcols())
    throw sgpp::base::algorithm_exception(
        "Row or column number of base matrix and result matrix are unequal.");

  // If level vectors are equal, no permutation has to be applied
  // Note: This leads to unnecessarily copying the unchanged base matrix and should be avoided.
  if (baseGridConfig.levelVector_ == desiredGridCOnfig.levelVector_) {
    permutatedMatrix = baseMatrix;
    return;
  }

  // Base level vector
  std::vector<size_t> baseLevelVec =
      PermutationUtil::deleteOnesFromLevelVec(baseGridConfig.levelVector_);
  std::vector<size_t> desiredLevelVec =
      PermutationUtil::deleteOnesFromLevelVec(desiredGridCOnfig.levelVector_);

  // Precomputation for h_inverse
  std::vector<size_t> preComputation = preComputeMatrixIndexForPoint(baseLevelVec);

  // compute permutation
  std::vector<size_t> permutation = computePermutation(baseLevelVec, desiredLevelVec);

  // TODO: Further validation

  // Permutation
  // First row is never permutated
  if (permutateRowsOrColums) {
    sgpp::base::DataVector firstRow(baseMatrix.getNcols());
    baseMatrix.getRow(0, firstRow);
    permutatedMatrix.setRow(0, firstRow);
  } else {
    sgpp::base::DataVector firstColumn(baseMatrix.getNrows());
    baseMatrix.getColumn(0, firstColumn);
    permutatedMatrix.setColumn(0, firstColumn);
  }
  int row = 1;
  // Init index and level vector
  std::vector<size_t> index(desiredLevelVec.size(), 1);
  std::vector<size_t> level(desiredLevelVec.size(), 1);

  // struct to reprensent a index/ level vector pair
  struct Point {
    std::vector<size_t> index;
    std::vector<size_t> level;
  };

  // already handled points
  std::vector<Point> points;
  // Add first point which is already finished
  Point firstPoint;
  firstPoint.index = index;
  firstPoint.level = level;
  points.push_back(firstPoint);

  for (size_t dim = 0; dim < desiredLevelVec.size(); dim++) {
    size_t currentSize = points.size();
    // Iterate over all current points
    for (int p = 0; p < currentSize; p++) {
      bool first = true;
      for (size_t l = 1; l <= desiredLevelVec[dim]; l++) {
        for (size_t i = 1; i < static_cast<size_t>(1 << l); i += 2) {
          // Element of level 1 has already been handled
          if (first) {
            first = false;
          } else {
            // Add new point to handled points
            Point newPoint;
            newPoint.level = points[p].level;
            newPoint.index = points[p].index;
            newPoint.level[dim] = l;
            newPoint.index[dim] = i;
            points.push_back(newPoint);
            // Get corresponding matrix row from old Q
            /*std::vector<size_t> baseLevel =
                PermutationUtil::permutateVector(newPoint.level, desiredLevelVec, baseLevelVec);
            std::vector<size_t> baseIndex =
                PermutationUtil::permutateVector(newPoint.index, desiredLevelVec, baseLevelVec);*/
            std::vector<size_t> baseLevel = applyPermutation(newPoint.level, permutation);
            std::vector<size_t> baseIndex = applyPermutation(newPoint.index, permutation);

            size_t correspondingBaseRowIndex =
                getMatrixIndexForPoint(baseLevel, baseIndex, baseLevelVec, preComputation) - 1;

            if (permutateRowsOrColums) {
              sgpp::base::DataVector baseRow(baseMatrix.getNcols());
              baseMatrix.getRow(correspondingBaseRowIndex, baseRow);
              permutatedMatrix.setRow(row, baseRow);
            } else {
              sgpp::base::DataVector baseColumn(baseMatrix.getNrows());
              baseMatrix.getColumn(correspondingBaseRowIndex, baseColumn);
              permutatedMatrix.setColumn(row, baseColumn);
            }
            // Increment row counter
            row++;
          }
        }
      }
    }
  }
}

void DBMatOfflinePermutable::dimensionBlowUp(
    const sgpp::base::GeneralGridConfiguration& baseGridConfig,
    const sgpp::base::GeneralGridConfiguration& desiredGridCOnfig,
    sgpp::base::DataMatrix& baseMatrix, bool matrixIsInverse) {
  auto& gridType = desiredGridCOnfig.type_;
  int dimDelta = desiredGridCOnfig.dim_ - baseGridConfig.dim_;
  // no dimension blow-up needs to be applied
  if (dimDelta == 0) return;
  if (gridType == sgpp::base::GridType::Linear) {
    double dimFactor = std::pow(3, std::abs(dimDelta));
    if (matrixIsInverse) {
      if (dimDelta > 0)
        baseMatrix.mult(dimFactor);
      else
        baseMatrix.mult(1.0 / dimFactor);
    } else {
      if (dimDelta > 0)
        baseMatrix.mult(1.0 / dimFactor);
      else
        baseMatrix.mult(dimFactor);
    }
  }
}

void DBMatOfflinePermutable::permuteLhsMatrix(
    const sgpp::base::GeneralGridConfiguration& baseGridConfig,
    const sgpp::base::GeneralGridConfiguration& desiredGridCOnfig) {
  // Copy base matrix for permutation
  sgpp::base::DataMatrix baseLhs(this->lhsMatrix);
  // Permutate rows
  permuteMatrix(baseGridConfig, desiredGridCOnfig, baseLhs, this->lhsMatrix, true);
  // Copy permutated lhs
  sgpp::base::DataMatrix rowPermutatedLhs(this->lhsMatrix);
  // Permutate collumns
  permuteMatrix(baseGridConfig, desiredGridCOnfig, rowPermutatedLhs, this->lhsMatrix, false);
  // Multiply dimensios blow up factor
  dimensionBlowUp(baseGridConfig, desiredGridCOnfig, this->lhsMatrix);
}

}  // namespace datadriven
}  // namespace sgpp