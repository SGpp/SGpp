/*
 * CombiMinMaxSchemeFT.hpp
 *
 *  Created on: Oct 19, 2015
 *      Author: sccs
 */

#ifndef SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEMEFT_HPP_
#define SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEMEFT_HPP_

#include "../../distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"

namespace combigrid {

class CombiMinMaxSchemeFT: public CombiMinMaxScheme {
public:

  CombiMinMaxSchemeFT(DimType dim, LevelVector& lmin, LevelVector& lmax) :
      CombiMinMaxScheme(dim, lmin, lmax) {

    // Add extra combiSpaces to ensure fault tolerance
    for (size_t i = 0; i < levels_.size(); ++i) {
      LevelVector& l = levels_[i];
      for (LevelType p = LevelType(effDim_);
          p < LevelType(effDim_) + extraDiags_; ++p) {
        if (l >= lmin && sum(l) == n_ - p) {
          combiSpaces_.push_back(l);
          coefficients_.push_back(0.0);
        }
      }
    }
  }

  virtual ~CombiMinMaxSchemeFT() {
  }

private:

  // Number of extra diagonals for fault tolerance. Recommended = 2
  const int extraDiags_ = 2;
};

} /* namespace combigrid */

#endif /* SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEMEFT_HPP_ */
