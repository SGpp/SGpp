/*
 * LevelVector.hpp
 *
 *  Created on: May 14, 2013
 *      Author: heenemo
 */

#ifndef INDEXVECTOR_HPP_
#define INDEXVECTOR_HPP_

#include <vector>
#include <assert.h>
#include <ostream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "sgpp/distributedcombigrid/utils/Types.hpp"

namespace combigrid {

typedef std::vector<IndexType> IndexVector;

inline IndexType sum(const IndexVector &l) {
  IndexType lsum(0);
  for (std::size_t d = 0; d < l.size(); ++d)
    lsum += l[d];

  return lsum;
}

inline bool operator==(const IndexVector &l1, const IndexVector &l2) {
  assert(l1.size() == l2.size());

  for (std::size_t i = 0; i < l1.size(); ++i) {
    if (l1[i] != l2[i])
      return false;
  }

  return true;
}

// a l1 < l2 if each entry l1,i < l2,i
inline bool operator<=(const IndexVector &l1, const IndexVector &l2) {
  assert(l1.size() == l2.size());

  for (std::size_t i = 0; i < l1.size(); ++i) {
    if (l1[i] > l2[i])
      return false;
  }

  return true;
}

// a l1 < l2 if each entry l1,i < l2,i
inline bool operator<(const IndexVector &l1, const IndexVector &l2) {
  assert(l1.size() == l2.size());

  for (std::size_t i = 0; i < l1.size(); ++i) {
    if (!(l1[i] < l2[i]))
      return false;
  }

  return true;
}

// a l1 > l2 if each entry l1,i > l2,i
inline bool operator>(const IndexVector &l1, const IndexVector &l2) {
  assert(l1.size() == l2.size());

  for (std::size_t i = 0; i < l1.size(); ++i) {
    if (!(l1[i] > l2[i]))
      return false;
  }

  return true;
}

// a l1 >= l2 if each entry l1,i >= l2,i
// todo replace some of these operators by using contrast operators
inline bool operator>=(const IndexVector &l1, const IndexVector &l2) {
  assert(l1.size() == l2.size());

  for (std::size_t i = 0; i < l1.size(); ++i) {
    if (l1[i] < l2[i])
      return false;
  }

  return true;
}

inline IndexVector operator+(const IndexVector &l1, const IndexVector &l2) {
  assert(l1.size() == l2.size());

  IndexVector tmp(l1.size());

  for (std::size_t i = 0; i < l1.size(); ++i)
    tmp[i] = l1[i] + l2[i];

  return tmp;
}

inline std::ostream& operator<<(std::ostream& os, const IndexVector& l) {
  os << "[";
  for (size_t i = 0; i < l.size(); ++i)
    os << l[i] << " ";
  os << "]";

  return os;
}

inline IndexVector operator-(const IndexVector &l1, const IndexVector &l2) {
  assert(l1.size() == l2.size());
  assert(l1 >= l2);

  IndexVector tmp(l1.size());

  for (std::size_t i = 0; i < l1.size(); ++i) {
    tmp[i] = l1[i] - l2[i];
  }

  return tmp;
}

/* read in indexvector from string where ' ' serves as delimiter
 * the vector described by str MUST NOT have a different size than ivec */
inline IndexVector& operator>>(std::string str, IndexVector& ivec) {
  std::vector<std::string> strs;
  boost::split(strs, str, boost::is_any_of(" "));

  assert(ivec.size() == strs.size());

  for (size_t i = 0; i < strs.size(); ++i)
    ivec[i] = boost::lexical_cast<int>(strs[i]);

  return ivec;
}

} // namespace combigrid

#endif /* LEVELVECTOR_HPP_ */
