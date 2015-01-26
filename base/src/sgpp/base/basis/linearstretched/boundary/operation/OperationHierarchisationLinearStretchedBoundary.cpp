/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de) */
/*        2010 - 2011 Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)    */
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


#include <sgpp/base/basis/linearstretched/boundary/operation/OperationHierarchisationLinearStretchedBoundary.hpp>
#include <sgpp/base/basis/linearstretched/boundary/algorithm_sweep/HierarchisationLinearStretchedBoundary.hpp>
#include <sgpp/base/basis/linearstretched/boundary/algorithm_sweep/DehierarchisationLinearStretchedBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>


namespace sg {
  namespace base {

    void OperationHierarchisationLinearStretchedBoundary::doHierarchisation(DataVector& node_values) {
      HierarchisationLinearStretchedBoundary func(this->storage);
      sweep<HierarchisationLinearStretchedBoundary> s(func, this->storage);

      // N D case
      if (this->storage->dim() > 1) {
        for (size_t i = 0; i < this->storage->dim(); i++) {
          s.sweep1D_Boundary(node_values, node_values, i);
        }
      }
      // 1 D case
      else {
        s.sweep1D(node_values, node_values, 0);
      }
    }

    void OperationHierarchisationLinearStretchedBoundary::doDehierarchisation(DataVector& alpha) {
      DehierarchisationLinearStretchedBoundary func(this->storage);
      sweep<DehierarchisationLinearStretchedBoundary> s(func, this->storage);

      // N D case
      if (this->storage->dim() > 1) {
        for (size_t i = 0; i < this->storage->dim(); i++) {
          s.sweep1D_Boundary(alpha, alpha, i);
        }
      }
      // 1 D case
      else {
        s.sweep1D(alpha, alpha, 0);
      }
    }

  }
}
