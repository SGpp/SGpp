#ifndef DMSYSTEMMATRIXSUBSPACES_HPP
#define DMSYSTEMMATRIXSUBSPACES_HPP

#include <string>

#include <base/datatypes/DataVector.hpp>
#include <base/grid/Grid.hpp>
#include <datadriven/algorithm/DMSystemMatrixBase.hpp>
#include "../operation/OperationMultipleEvalSubspace/CommonParameters.hpp"

#include "base/operation/OperationMultipleEval.hpp"

#include "datadriven/operation/DatadrivenOperationCommon.hpp"
#include "datadriven/DatadrivenOpFactory.hpp"

//#include "AbstractOperationMultipleEval.hpp"

namespace sg {
namespace datadriven {

/**
 * Class that implements the virtual class sg::base::OperationMatrix for the
 * application of classification for the Systemmatrix
 *
 * The Identity matrix is used as regularization operator.
 *
 * For the Operation B's mult and mutlTransposed functions
 * vectorized formulations are used.
 */
class SystemMatrixLeastSquaresIdentity: public sg::datadriven::DMSystemMatrixBase {
private:
	/// vectorization mode
	//ComputeKernelType kernelType;
	/// Number of orignal training instances
	size_t instances;
	/// Number of patched and used training instances
	size_t paddedInstances;
	/// OperationB for calculating the data matrix
	//AbstractOperationMultipleEval* B;
	sg::base::OperationMultipleEval *B;

	sg::base::Grid &grid;

	sg::datadriven::OperationMultipleEvalConfiguration implementationConfiguration;

public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param trainData reference to sg::base::DataMatrix that contains the training data
	 * @param lambda the lambda, the regression parameter
	 * @param kernelType compute kernel used
	 */
	SystemMatrixLeastSquaresIdentity(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda);

	/**
	 * Std-Destructor
	 */
	virtual ~SystemMatrixLeastSquaresIdentity();

	virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

	virtual void generateb(sg::base::DataVector& classes, sg::base::DataVector& b);

	virtual void rebuildLevelAndIndex();

	void setImplementation(sg::datadriven::OperationMultipleEvalConfiguration operationConfiguration) {
		this->implementationConfiguration = operationConfiguration;
		this->B = sg::op_factory::createOperationMultipleEval(this->grid, *(this->dataset_), this->implementationConfiguration);
	}
};

}
}

#endif /* DMSYSTEMMATRIXSUBSPACES_HPP */
