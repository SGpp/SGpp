#ifndef LEARNERDENSITYBASED_HPP_
#define LEARNERDENSITYBASED_HPP_

#include "datadriven/application/LearnerBase.hpp"
#include "datadriven/application/Learner.hpp"

namespace sg
{

namespace datadriven
{


class LearnerDensityBased: public sg::datadriven::LearnerBase
{
protected:
	//Mapping from class index to class number:
	std::map<int, double> index_to_class_;
	//Stores the coefficients for every class
	std::vector<sg::base::DataVector> alphas_;
	/// regularization operator
	sg::base::OperationMatrix* C_;
	/// regularization mode
	sg::datadriven::LearnerRegularizationType CMode_;
public:
	LearnerDensityBased(sg::datadriven::LearnerRegularizationType&, const bool isRegression, const bool isVerbose = true);
	virtual ~LearnerDensityBased();
	/**
	 * Learning a dataset with spatially adaptive sparse grids
	 *
	 * @param testDataset the training dataset
	 * @param classes classes corresponding to the training dataset
	 * @param GridConfig configuration of the regular start grid
	 * @param SolverConfigRefine configuration of the SLE solver during the adaptive refinements of the grid
	 * @param SolverConfigFinal configuration of the final SLE solving step on the refined grid
	 * @param AdaptConfig configuration of the adaptivity strategy
	 * @param testAccDuringAdapt set to true if the training accuracy should be determined in evert refinement step
	 * @param lambda regularization parameter lambda
	 */
	virtual LearnerTiming train(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes,
			const sg::base::RegularGridConfiguration& GridConfig, const sg::solver::SLESolverConfiguration& SolverConfigRefine,
			const sg::solver::SLESolverConfiguration& SolverConfigFinal, const sg::base::AdpativityConfiguration& AdaptConfig,
			bool testAccDuringAdapt, const double lambda);
	virtual sg::base::DataVector predict(sg::base::DataMatrix& testDataset);
	/// construct system matrix
		virtual sg::datadriven::DMSystemMatrixBase* createDMSystem(sg::base::DataMatrix& trainDataset, double lambda);
};

}

}

#endif /* LEARNERDENSITYBASED_HPP_ */
