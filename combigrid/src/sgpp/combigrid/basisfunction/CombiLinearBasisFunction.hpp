// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef COMBILINEARBASISFUNCTION_HPP_
#define COMBILINEARBASISFUNCTION_HPP_


#include <sgpp/combigrid/utils/combigrid_ultils.hpp>
#include <sgpp/combigrid/basisfunction/CombiBasisFunctionBasis.hpp>

namespace combigrid {

  /** Linear basis function */
  class LinearBasisFunction : public combigrid::BasisFunctionBasis {
    public:

      /** empty Ctror */
      LinearBasisFunction() {
        ;
      }

      /** first method which returns the contribution of the first point in the 1D cell
       * @param coord  1D coordonate idealy should be [0,1] but for extrapolation could be different [-1,2]*/
      virtual double functionEval1(double coord) const {
        return (1.0 - coord);
      }

      /** second method which returns the contribution of the second point in the 1D cell
       * @param coord  1D coordonate idealy should be [0,1] but for extrapolation could be different [-1,2]*/
      virtual double functionEval2(double coord) const {
        return (coord);
      }

      /** return the default basis function*/
      static const BasisFunctionBasis* getDefaultBasis() {
        return defaultBasis_;
      }

    private:

      /** default basis function */
      static const BasisFunctionBasis* defaultBasis_;

  };
}

#endif /* COMBILINEARBASISFUNCTION_HPP_ */