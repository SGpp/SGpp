

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflinePermutable.hpp>

namespace sgpp {
namespace datadriven {

const std::set<MatrixDecompositionType> DBMatOfflinePermutable::PermutableDecompositions{
    MatrixDecompositionType::OrthoAdapt
    };

DBMatOfflinePermutable::DBMatOfflinePermutable() : DBMatOffline() {}

DBMatOfflinePermutable::DBMatOfflinePermutable(const std::string& fileName)
    : DBMatOffline(fileName) {}

size_t DBMatOfflinePermutable::getMatrixIndexForPoint(std::vector<size_t> level,
                                                      std::vector<size_t> index,
                                                      std::vector<size_t> gridLevel) {
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

    size_t prod = 1;
    for (size_t i = 0; i < lStar; i++) {
      prod *= (1 << gridLevel[i]) - 1;
    }

    result += (prod + (1 << (level[lStar] - 1)) + ((index[lStar] + 1) >> 1) - 3) * mult;

    mult *= (1 << gridLevel[lStar]) - 2;

    // determine new lStar
    int lStar_ = -1;
    for (size_t i = lStar - 1; i >=  0; i--) {
      if (level[i] > 1) lStar_ = i;
    }
    lStar = lStar_;
  }
  return result;
}

std::vector<size_t> DBMatOfflinePermutable::deleteOnesFromLevelVec(
    std::vector<size_t> vectorWithOnes) {
  std::vector<size_t> output;
  for (size_t i = 0; i < vectorWithOnes.size(); i++) {
    if (vectorWithOnes[i] != 1) output.push_back(vectorWithOnes[i]);
  }
  return output;
}

std::vector<size_t> DBMatOfflinePermutable::permutateVector(std::vector<size_t> vector,
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

void DBMatOfflinePermutable::permutateMatrix(sgpp::base::CombiGridConfiguration baseGridConfig,
                                             sgpp::base::CombiGridConfiguration desiredGridCOnfig,
                                             const sgpp::base::DataMatrix& baseMatrix,
                                             sgpp::base::DataMatrix& permutatedMatrix,
                                             bool permutateRowsOrColums) {
  // Base level vector
  std::vector<size_t> baseLevelVec = baseGridConfig.levels;

  // TODO: Further validation

  std::vector<size_t> desiredLevelVec = deleteOnesFromLevelVec(desiredGridCOnfig.levels);

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
            std::vector<size_t> baseLevel =
                permutateVector(newPoint.level, desiredLevelVec, baseLevelVec);
            std::vector<size_t> baseIndex =
                permutateVector(newPoint.index, desiredLevelVec, baseLevelVec);

            size_t correspondingBaseRowIndex =
                getMatrixIndexForPoint(baseLevel, baseIndex, baseLevelVec) - 1;

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

void DBMatOfflinePermutable::dimensionBlowUp(sgpp::base::CombiGridConfiguration baseGridConfig,
                                             sgpp::base::CombiGridConfiguration desiredGridCOnfig,
                                             sgpp::base::DataMatrix& baseMatrix,
                                             bool matrixIsInverse) {
  std::cout << "dimension blow up..."
            << "\n";
  auto& gridType = desiredGridCOnfig.type_;
  int dimDelta = desiredGridCOnfig.dim_ - baseGridConfig.dim_;
  if (gridType == sgpp::base::GridType::Linear) {
    float dimFactor = std::pow(3, std::abs(dimDelta));
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

void DBMatOfflinePermutable::permutateLhsMatrix(
    sgpp::base::CombiGridConfiguration baseGridConfig,
    sgpp::base::CombiGridConfiguration desiredGridCOnfig) {
  // Copy base matrix for permutation
  sgpp::base::DataMatrix baseLhs(this->lhsMatrix);
  // Permutate rows
  permutateMatrix(baseGridConfig, desiredGridCOnfig, baseLhs, this->lhsMatrix, true);
  // Copy permutated lhs
  sgpp::base::DataMatrix rowPermutatedLhs(this->lhsMatrix);
  // Permutate collumns
  permutateMatrix(baseGridConfig, desiredGridCOnfig, rowPermutatedLhs, this->lhsMatrix, false);
  // Multiply dimensios blow up factor
  dimensionBlowUp(baseGridConfig, desiredGridCOnfig, this->lhsMatrix);
}

}  // namespace datadriven
}  // namespace sgpp