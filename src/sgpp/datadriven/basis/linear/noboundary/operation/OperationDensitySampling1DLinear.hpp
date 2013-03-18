/*
 * OperationDensitySampling1DLinear.hpp
 *
 *  Created on: Dec 4, 2012
 *      Author: Emily Mo-Hellenbrand
 */

#ifndef OPERATIONDENSITYSAMPLING1DLINEAR_HPP_
#define OPERATIONDENSITYSAMPLING1DLINEAR_HPP_

#include "base/grid/Grid.hpp"
#include "datadriven/operation/OperationDensitySampling1D.hpp"

namespace sg
{
namespace datadriven
{
	class OperationDensitySampling1DLinear : public OperationDensitySampling1D
	{
		protected:
			base::Grid* grid;
		public:
			OperationDensitySampling1DLinear(base::Grid* grid);
			virtual ~OperationDensitySampling1DLinear();

		    /**
		     * Sampling on 1D grid
		     *
		     * @param alpha Coefficient vector for current grid (1D grid)
		     * @param num_samples # of samples to draw
		     * @param samples Output DataVector
		     */
			void doSampling1D(base::DataVector* alpha, size_t num_samples, base::DataVector* &samples);
	};

}
}

#endif /* OPERATIONDENSITYSAMPLING1D_HPP_ */
