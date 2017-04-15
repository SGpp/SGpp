/*
 * DBMatOfflineLU.hpp
 *
 *  Created on: 02.03.2017
 *      Author: michael
 */

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOfflineGE.hpp>

namespace sgpp {
namespace datadriven {

class DBMatOfflineLU : public DBMatOfflineGE {
 public:
  DBMatOfflineLU(const DBMatDensityConfiguration& oc);

  DBMatOfflineLU(const DBMatOfflineLU& rhs);

  DBMatOfflineLU(DBMatOfflineLU&& rhs) = default;

  DBMatOfflineLU& operator=(const DBMatOfflineLU& rhs);

  DBMatOfflineLU& operator=(DBMatOfflineLU&& rhs) = default;

  virtual ~DBMatOfflineLU() = default;

  DBMatOffline* clone() override;

  bool isRefineable() override;

  DBMatOfflineLU(const std::string& fileName);

  void decomposeMatrix() override;

  void permuteVector(DataVector& b);

  /**
   * Store the decomposed matrix, the permutation and configuration.
   *
   * @param fname the file name
   */
  void store(const std::string& fname) override;

 private:
  std::unique_ptr<gsl_permutation> permutation;  // Stores the permutation that was
                                                 // applied on the matrix during decomposition
};

} /* namespace datadriven */
} /* namespace sgpp */
