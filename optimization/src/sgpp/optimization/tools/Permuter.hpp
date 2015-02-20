// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TOOLS_PERMUTER_HPP
#define SGPP_OPTIMIZATION_TOOLS_PERMUTER_HPP

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Helper class for sorting vectors.
     *
     * This class is designed to make the permutation of a vector \f$\vec{a}\f$
     * according to another vector \f$\vec{b}\f$ possible.
     *
     * \f$\vec{a} = (a_1, \dotsc, a_n)\f$ can be reordered to
     * \f$(a_{\pi(1)}, \dotsc, a_{\pi(n)})\f$ with
     * \f$\pi\f$ a permutation of \f$\{1, \dotsc, n\}\f$ such that
     * \f$b_{\pi(1)} < \dotsb < b_{\pi(n)}\f$
     * (realising permutation of the sorting of \f$\vec{b}\f$),
     * assuming \f$\vec{b}\f$ contains no equal elements.
     *
     * The type T must have a "<" operator.
     *
     * Example usage:
     *
     * \code{.cpp}
     * #include <algorithm>
     * std::vector<float_t> a, b;
     * ...
     * Permuter<float_t> permuter(b);
     * std::sort(a.begin(), a.end(), permuter);
     * \endcode
     */
    template <class T>
    class Permuter {
      public:
        /**
         * Constructor.
         * Don't destruct vec before this object!
         *
         * @param vec   reference to \f$\vec{b}\f$
         */
        Permuter(const std::vector<T>& vec) : vec(vec) {
        }

        /**
         * Comparison operator for std::sort.
         *
         * @param a     index 1
         * @param b     index 2
         * @return      (vec[a] < vec[b])
         */
        inline bool operator()(const size_t& a, const size_t& b) const {
          return (vec[a] < vec[b]);
        }

      protected:
        /// reference to \f$\vec{b}\f$
        const std::vector<T>& vec;
    };

  }
}

#endif /* SGPP_OPTIMIZATION_TOOLS_PERMUTER_HPP */
