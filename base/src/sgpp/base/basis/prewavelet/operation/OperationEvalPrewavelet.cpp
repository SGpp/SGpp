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
#include <sgpp/base/basis/prewavelet/operation/OperationEvalPrewavelet.hpp>

#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>



#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    double OperationEvalPrewavelet::eval(DataVector& alpha,
                                         std::vector<double>& point) {
      typedef std::vector<std::pair<size_t, double> > IndexValVector;

      IndexValVector vec;
      PrewaveletBasis<unsigned int, unsigned int> base;
      GetAffectedBasisFunctions<PrewaveletBasis<unsigned int, unsigned int> > ga(
        storage);

      ga(base, point, vec);

      double result = 0.0;

      for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
        result += iter->second * alpha[iter->first];
      }

      return result;
    }

    double OperationEvalPrewavelet::test(DataVector& alpha, DataVector& data,
                                         DataVector& classes) {
      return 0;
    }

    double OperationEvalPrewavelet::integrate(DataVector& alpha) {
      double result = 0.0;

      for (size_t i = 0; i < storage->size(); i++) {
        double temp_result = 1;

        for (size_t d = 0; d < storage->dim(); d++) {
          GridStorage::index_type::level_type level;
          GridStorage::index_type::index_type index;
          (*storage)[i]->get(d, level, index);

          if (index != 1 && index != (unsigned int)((1 << level) - 1)) {
            temp_result = 0.0;
            break;
          } else if (level == 1) {
            temp_result *= 1.0 / 2.0;
          } else {
            temp_result *= 0.4 * (1.0 / (1 << level));
          }

        }

        result += alpha[i] * temp_result;
      }

      return result;
    }

  }
}
