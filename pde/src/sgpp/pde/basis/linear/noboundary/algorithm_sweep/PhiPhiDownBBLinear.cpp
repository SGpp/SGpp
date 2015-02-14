// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {



    PhiPhiDownBBLinear::PhiPhiDownBBLinear(SGPP::base::GridStorage* storage) : storage(storage), boundingBox(storage->getBoundingBox()) {
    }

    PhiPhiDownBBLinear::~PhiPhiDownBBLinear() {
    }

    void PhiPhiDownBBLinear::operator()(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim) {
      //std::cout << dim << std::endl;
      //std::cout << index.toString() << std::endl;

      float_t q = this->boundingBox->getIntervalWidth(dim);
      float_t t = this->boundingBox->getIntervalOffset(dim);

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

    void PhiPhiDownBBLinear::rec(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim, float_t fl, float_t fr) {
      size_t seq = index.seq();

      float_t alpha_value = source[seq];

      SGPP::base::GridStorage::index_type::level_type l;
      SGPP::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      float_t h = 1.0 / static_cast<float_t>(1 << l);
      float_t tmp_m = ((fl + fr) / 2.0);

      // integration
      result[seq] = (h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value);

      // dehierarchisation
      float_t fm = tmp_m + alpha_value;

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

    void PhiPhiDownBBLinear::recBB(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim, float_t fl, float_t fr, float_t q, float_t t) {
      size_t seq = index.seq();

      float_t alpha_value = source[seq];

      SGPP::base::GridStorage::index_type::level_type l;
      SGPP::base::GridStorage::index_type::index_type i;

      index.get(dim, l, i);

      float_t h = 1.0 / static_cast<float_t>(1 << l);
      float_t tmp_m = ((fl + fr) / 2.0);

      // integration
      result[seq] = ((h * tmp_m) + (((2.0 / 3.0) * h) * alpha_value)) * q;  // diagonal entry

      // dehierarchisation
      float_t fm = tmp_m + alpha_value;

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
  }
  // namespace SGPP
}