/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
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

#include <sgpp/base/basis/prewavelet/algorithm_sweep/ConvertLinearToPrewavelet.hpp>
#include <sgpp/base/basis/prewavelet/algorithm_sweep/ConvertPrewaveletToLinear.hpp>
#include <sgpp/base/basis/prewavelet/operation/OperationConvertPrewavelet.hpp>


#include <sgpp/base/algorithm/sweep.hpp>



namespace sg {
  namespace base {

    void OperationConvertPrewavelet::doConvertToLinear(
      DataVector& alpha) {

      ConvertPrewaveletToLinear func(this->storage);
      sweep<ConvertPrewaveletToLinear> s(func, this->storage);


      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(alpha, alpha, i);
      }

    }

    void OperationConvertPrewavelet::doConvertFromLinear(DataVector& alpha) {
      ConvertLinearToPrewavelet func(this->storage, this->shadowstorage);
      sweep<ConvertLinearToPrewavelet> s(func, this->storage);

      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(alpha, alpha, i);
      }
    }

  }
}
