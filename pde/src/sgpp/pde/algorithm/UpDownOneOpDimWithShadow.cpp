/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#include <sgpp/pde/algorithm/UpDownOneOpDimWithShadow.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

namespace sg {
  namespace pde {

    UpDownOneOpDimWithShadow::UpDownOneOpDimWithShadow(sg::base::GridStorage* storage,
        sg::base::GridStorage* shadowStorage) {
      this->storage = storage;
      this->shadowStorage = shadowStorage;
    }

    UpDownOneOpDimWithShadow::~UpDownOneOpDimWithShadow() {
    }

    void UpDownOneOpDimWithShadow::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {

      expandGrid();

      //Create new Datavectors for the grid including shadow points.
      sg::base::DataVector alpha_temp(storage->size());
      sg::base::DataVector result_temp(storage->size());

      alpha_temp.setAll(0.0);
      result_temp.setAll(0.0);

      //Copy the alpha vector ... the old points remain the same, the new points are zero.
      for (size_t i = 0; i < alpha.getSize(); i++) {
        alpha_temp[i] = alpha[i];
      }

      sg::base::DataVector beta(result_temp.getSize());
      beta.setAll(0.0);

      for (size_t i = 0; i < storage->dim(); i++) {
        this->updown(alpha_temp, beta, storage->dim() - 1, i);

        result_temp.add(beta);

      }

      //Remove shadow grid points from the grid
      shrinkGrid();

      //Copy the result in the actual result vector. The values for the shadow
      //grid points are just needed for data transport, thus they can be omitted
      for (size_t i = 0; i < storage->size(); i++) {
        result[i] = result_temp[i];
      }
    }

    void UpDownOneOpDimWithShadow::expandGrid() {
      for (size_t i = 0; i < shadowStorage->size(); i++) {
        storage->insert(*shadowStorage->get(i));
      }
    }

    void UpDownOneOpDimWithShadow::shrinkGrid() {
      for (size_t i = 0; i < shadowStorage->size(); i++) {
        storage->deleteLast();
      }
    }

    void UpDownOneOpDimWithShadow::updown(sg::base::DataVector& alpha, sg::base::DataVector& result,
                                          size_t dim, size_t op_dim) {
      if (dim == op_dim) {
        specialOP(alpha, result, dim, op_dim);
      } else {
        //Unidirectional scheme
        if (dim > 0) {
          // Reordering ups and downs
          sg::base::DataVector temp(alpha.getSize());
          temp.setAll(0.0);
          up(alpha, temp, dim);
          updown(temp, result, dim - 1, op_dim);

          // Same from the other direction:
          sg::base::DataVector result_temp(alpha.getSize());
          result_temp.setAll(0.0);
          updown(alpha, temp, dim - 1, op_dim);
          down(temp, result_temp, dim);

          result.add(result_temp);
        } else {
          // Terminates dimension recursion
          up(alpha, result, dim);

          sg::base::DataVector temp(alpha.getSize());
          temp.setAll(0.0);
          down(alpha, temp, dim);

          result.add(temp);
        }

      }
    }

    void UpDownOneOpDimWithShadow::specialOP(sg::base::DataVector& alpha, sg::base::DataVector& result,
        size_t dim, size_t op_dim) {
      //Unidirectional scheme
      if (dim > 0) {
        // Reordering ups and downs
        sg::base::DataVector temp(alpha.getSize());
        temp.setAll(0.0);
        upOpDim(alpha, temp, dim);
        updown(temp, result, dim - 1, op_dim);

        // Same from the other direction:
        sg::base::DataVector result_temp(alpha.getSize());
        result_temp.setAll(0.0);
        updown(alpha, temp, dim - 1, op_dim);
        downOpDim(temp, result_temp, dim);

        result.add(result_temp);
      } else {
        // Terminates dimension recursion
        upOpDim(alpha, result, dim);

        sg::base::DataVector temp(alpha.getSize());
        temp.setAll(0.0);
        downOpDim(alpha, temp, dim);

        result.add(temp);
      }
    }

  }
}
