/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#include <sgpp/base/basis/modpoly/ModifiedPolyBasis.hpp>

#include <sgpp/base/basis/prewavelet/operation/OperationMultipleEvalPrewavelet.hpp>

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    void OperationMultipleEvalPrewavelet::mult(DataVector& alpha, DataVector& result) {
      AlgorithmDGEMV<SPrewaveletBase> op;
      PrewaveletBasis<unsigned int, unsigned int> base;

      op.mult(storage, base, alpha, this->dataset, result);
    }

    void OperationMultipleEvalPrewavelet::multTranspose(DataVector& source, DataVector& result) {
      AlgorithmDGEMV<SPrewaveletBase> op;
      PrewaveletBasis<unsigned int, unsigned int> base;

      op.mult_transposed(storage, base, source, this->dataset, result);
    }

  }
}
