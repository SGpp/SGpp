/*
 * CombiLinearBasisFunction.hpp
 *
 *  Created on: Feb 21, 2011
 *      Author: benk
 */

#ifndef COMBILINEARBASISFUNCTION_HPP_
#define COMBILINEARBASISFUNCTION_HPP_


#include "combigrid/utils/combigrid_ultils.hpp"
#include "combigrid/basisfunction/CombiBasisFunctionBasis.hpp"

namespace combigrid{

   /** Linear basis function */
   class LinearBasisFunction : public combigrid::BasisFunctionBasis {
   public:

	   /** empty Ctror */
	   LinearBasisFunction() {;}

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
	   static const BasisFunctionBasis* getDefaultBasis() { return defaultBasis_;}

   private:

	   /** default basis function */
	   static const BasisFunctionBasis* defaultBasis_;

   };
}

#endif /* COMBILINEARBASISFUNCTION_HPP_ */
