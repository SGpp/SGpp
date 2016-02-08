/*
 * CombiMinMaxScheme.hpp
 *
 *  Created on: Oct 2, 2015
 *      Author: sccs
 */

#ifndef SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_
#define SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_

#include <boost/math/special_functions/binomial.hpp>
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"

namespace combigrid {

class CombiMinMaxScheme {
public:
  CombiMinMaxScheme(DimType dim, LevelVector& lmin, LevelVector& lmax) {
    assert(dim > 0);

    assert(lmax.size() == dim);
    for (size_t i = 0; i < lmax.size(); ++i)
      assert(lmax[i] > 0);

    assert(lmin.size() == dim);
    for (size_t i = 0; i < lmin.size(); ++i)
      assert(lmin[i] > 0);

    dim_ = dim;
    // create hierarchical subspaces
    createLevels(dim, lmax, lmin);
    // create combi spaces
    for (size_t i = 0; i < levels_.size(); ++i) {
      LevelVector& l = levels_[i];
      for (LevelType p = 0; p < LevelType(effDim_); ++p) {
        if (l >= lmin && sum(l) == n_ - p) {
          combiSpaces_.push_back(l);
          coefficients_.push_back(
              std::pow(-1, p)
                  * boost::math::binomial_coefficient<double>(effDim_ - 1, p));
        }
      }
    }
  }

  ~CombiMinMaxScheme() {
  }

  inline const std::vector<LevelVector>& getCombiSpaces() const {
    return combiSpaces_;
  }

  inline const std::vector<double>& getCoeffs() const {
    return coefficients_;
  }

  inline void print(std::ostream& os) const;

protected:
  LevelType n_;

  DimType dim_;

  DimType effDim_;

  std::vector<LevelVector> levels_;

  std::vector<LevelVector> combiSpaces_;

  std::vector<double> coefficients_;

  void createLevels(DimType dim, const LevelVector& nmax, LevelVector& lmin);
  void createLevelsRec(DimType dim, LevelType n, DimType d, LevelVector &l,
      const LevelVector& nmax);
};

void CombiMinMaxScheme::createLevels(DimType dim, const LevelVector& nmax,
    LevelVector& lmin) {
  assert(nmax.size() == dim);

  // Remove dummy dimensions (e.g. if lmin = (2,2,2) and nmax = (4,4,2), dimension 3 is dummy)
  LevelVector nmaxtmp = nmax;
  LevelVector lmintmp = lmin;
  LevelVector dummyDims;
  for (size_t i = 0; i < nmax.size(); ++i) {
    if (nmax[i] - lmin[i] == 0) {
      nmaxtmp.erase(nmaxtmp.begin() + i - dummyDims.size());
      lmintmp.erase(lmintmp.begin() + i - dummyDims.size());
      dummyDims.push_back(i);
    }
  }

  effDim_ = nmaxtmp.size();
  // compute c which fulfills nmax - c*1  >= lmin
  LevelType c(0);
  LevelVector cv = nmaxtmp - lmintmp;

  // check if all elements in nmax-lmin are equal. If not, define c as the max of cv
  if (std::adjacent_find(cv.begin(), cv.end(), std::not_equal_to<int>())
      == cv.end()) {
    c = cv[effDim_ - 1];
  } else {
    c = *std::min_element(cv.begin(), cv.end());
  }
  LevelVector rlmin(effDim_);
  for (size_t i = 0; i < rlmin.size(); ++i) {
    rlmin[i] = nmaxtmp[i] - c;
  }

  // Define new lmin (dummy dimensions will be fixed later)
  for (size_t i = 0; i < lmin.size(); ++i) {
    lmin[i] = nmax[i] - c;
  }
  LevelType n = sum(rlmin) + c;

  LevelVector l(effDim_);
  createLevelsRec(effDim_, n, effDim_, l, nmaxtmp);

  // re-insert dummy dimensions left out
  for (size_t i = 0; i < dummyDims.size(); ++i) {
    for (size_t j = 0; j < levels_.size(); ++j) {
      levels_[j].insert(levels_[j].begin() + dummyDims[i], nmax[i]);
    }
    lmin[dummyDims[i]] = nmax[i];
  }
  n_ = sum(lmin) + c;
}

void CombiMinMaxScheme::createLevelsRec(DimType dim, LevelType n, DimType d,
    LevelVector &l, const LevelVector& nmax) {
  // sum rightmost entries of level vector
  LevelType lsum(0);
  for (size_t i = dim; i < l.size(); ++i)
    lsum += l[i];

  for (LevelType ldim = 1; ldim <= LevelType(n) + LevelType(d) - 1 - lsum;
      ++ldim) {
    l[dim - 1] = ldim;
    if (dim == 1) {
      if (l <= nmax) {
        levels_.push_back(l);
      }
    } else {
      createLevelsRec(dim - 1, n, d, l, nmax);
    }
  }
}

void CombiMinMaxScheme::print(std::ostream& os) const {
  for (uint i = 0; i < combiSpaces_.size(); ++i)
    os << "\t" << combiSpaces_[i] << "\t" << coefficients_[i] << std::endl;

  os << std::endl;
}

inline std::ostream& operator<<(std::ostream& os,
    const combigrid::CombiMinMaxScheme& scheme) {
  scheme.print(os);
  return os;
}

}
#endif /* SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_ */
