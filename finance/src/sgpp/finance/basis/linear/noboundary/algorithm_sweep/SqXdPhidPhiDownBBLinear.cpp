// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/SqXdPhidPhiDownBBLinear.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {



    SqXdPhidPhiDownBBLinear::SqXdPhidPhiDownBBLinear(SGPP::base::GridStorage* storage) : storage(storage), boundingBox(storage->getBoundingBox()) {
    }

    SqXdPhidPhiDownBBLinear::~SqXdPhidPhiDownBBLinear() {
    }

    void SqXdPhidPhiDownBBLinear::operator()(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim) {
      float_t q = boundingBox->getIntervalWidth(dim);
      float_t t = boundingBox->getIntervalOffset(dim);

      bool useBB = false;

      if (q != 1.0 || t != 0.0) {
        useBB = true;
      }

      if (useBB) {
        recBB(source, result, index, dim, 0.0, 0.0, q, t);
      } else {
        rec(source, result, index, dim, 0.0, 0.0);
      }
    }

    void SqXdPhidPhiDownBBLinear::rec(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim, float_t fl, float_t fr) {
      size_t seq = index.seq();

      float_t alpha_value = source[seq];

      SGPP::base::GridStorage::index_type::level_type l;
      SGPP::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      float_t i_dbl = static_cast<float_t>(i);
      int l_int = static_cast<int>(l);

      float_t diagonal = ((1.0 / 3.0) + (i_dbl * i_dbl)) * pow(2.0, 1 - l_int);

      // integration
      result[seq] = (  (((1.0 / static_cast<float_t>(1 << l_int)) * i_dbl) * (fl - fr)) + (diagonal * alpha_value) );

      // dehierarchisation
      float_t fm = ((fl + fr) / 2.0) + alpha_value;

      if (!index.hint()) {
        index.left_child(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, fl, fm);
        }

        index.step_right(dim);

        if (!storage->end(index.seq())) {
          rec(source, result, index, dim, fm, fr);
        }

        index.up(dim);
      }
    }

    void SqXdPhidPhiDownBBLinear::recBB(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim, float_t fl, float_t fr, float_t q, float_t t) {
      size_t seq = index.seq();

      float_t alpha_value = source[seq];

      SGPP::base::GridStorage::index_type::level_type l;
      SGPP::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      float_t i_dbl = static_cast<float_t>(i);
      int l_int = static_cast<int>(l);

      float_t diagonal = (1.0 / 3.0) * ((((pow(2.0, (1 - l_int))) * q * q) * (3.0 * (i_dbl * i_dbl) + 1)) + (12.0 * t * q * i_dbl) + (3.0 * (static_cast<float_t>(1 << (1 + l_int))) * t * t)) / (q);

      // integration
      result[seq] = (  (((1.0 / static_cast<float_t>(1 << l_int)) * i_dbl * q + t) * (fl - fr)) + (diagonal * alpha_value) );

      // dehierarchisation
      float_t fm = ((fl + fr) / 2.0) + alpha_value;

      if (!index.hint()) {
        index.left_child(dim);

        if (!storage->end(index.seq())) {
          recBB(source, result, index, dim, fl, fm, q, t);
        }

        index.step_right(dim);

        if (!storage->end(index.seq())) {
          recBB(source, result, index, dim, fm, fr, q, t);
        }

        index.up(dim);
      }
    }

    // namespace detail

  } // namespace SGPP
}