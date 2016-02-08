/*
 * CombiParameters.hpp
 *
 *  Created on: Dec 8, 2015
 *      Author: heenemo
 */

#ifndef SRC_SGPP_COMBIGRID_MANAGER_COMBIPARAMETERS_HPP_
#define SRC_SGPP_COMBIGRID_MANAGER_COMBIPARAMETERS_HPP_

#include <boost/serialization/map.hpp>
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"

namespace combigrid {

class CombiParameters {
 public:
  CombiParameters() {
  }

  CombiParameters(DimType dim, LevelVector lmin, LevelVector lmax,
                  std::vector<bool>& boundary, std::vector<LevelVector>& levels,
                  std::vector<real>& coeffs, std::vector<int>& taskIDs ) :
    dim_(dim), lmin_(lmin), lmax_(lmax), boundary_(boundary) {
    setLevelsCoeffs( taskIDs, levels, coeffs );
  }

  ~CombiParameters() {
  }

  inline const LevelVector& getLMin() {
    return lmin_;
  }

  inline const LevelVector& getLMax() {
    return lmax_;
  }

  inline const std::vector<bool>& getBoundary() {
    return boundary_;
  }

  inline real getCoeff(int taskID) {
    return coeffs_[taskID];
  }

  inline void getCoeffs(std::vector<int>& taskIDs, std::vector<real>& coeffs) {
    for (auto it : coeffs_) {
      taskIDs.push_back( it.first );
      coeffs.push_back( it.second );
    }
  }

  inline void setCoeff(int taskID, real coeff) {
    coeffs_[taskID] = coeff;
  }

  inline void setLevelsCoeffs(std::vector<int>& taskIDs,
                              std::vector<LevelVector>& levels, std::vector<real>& coeffs) {
    assert(taskIDs.size() == coeffs.size());
    assert(taskIDs.size() == levels.size());

    for (size_t i = 0; i < taskIDs.size(); ++i) {
      coeffs_[taskIDs[i]] = coeffs[i];
      levels_[taskIDs[i]] = levels[i];
    }
  }

  inline const LevelVector& getLevel( int taskID ) {
    return levels_[ taskID ];
  }

  inline DimType getDim() {
    return dim_;
  }

 private:
  DimType dim_;

  LevelVector lmin_;

  LevelVector lmax_;

  std::vector<bool> boundary_;

  std::map<int, LevelVector> levels_;

  std::map<int, real> coeffs_;

  friend class boost::serialization::access;

  // serialize
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version);
};

template<class Archive>
void CombiParameters::serialize(Archive& ar, const unsigned int version) {
  ar& dim_;
  ar& lmin_;
  ar& lmax_;
  ar& boundary_;
  ar& levels_;
  ar& coeffs_;
}

}

#endif /* SRC_SGPP_COMBIGRID_MANAGER_COMBIPARAMETERS_HPP_ */
