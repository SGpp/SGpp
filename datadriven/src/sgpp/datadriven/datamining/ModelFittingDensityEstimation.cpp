// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/ModelFittingDensityEstimation.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/exception/application_exception.hpp>

using namespace SGPP::base;

namespace SGPP {
namespace datadriven {

ModelFittingDensityEstimation::ModelFittingDensityEstimation(
		datadriven::DataMiningConfigurationDensityEstimation config) :
		datadriven::ModelFittingBase(), configuration(config) {
}

ModelFittingDensityEstimation::~ModelFittingDensityEstimation() {
}

void ModelFittingDensityEstimation::fit(datadriven::Dataset& dataset) {
	DataMatrix samples = dataset.getTrainingData();
	size_t numSamples = samples.getNrows();
	size_t numDims = samples.getNcols();

	initializeGrid(configuration.gridConfig);

	GridStorage* gridStorage = grid->getStorage();
	GridGenerator* gridGen = grid->createGridGenerator();
	DataVector rhs(grid->getStorage()->size());
	alpha->resize(grid->getStorage()->size());
	alpha->setAll(0.0);

	if (!configuration.sgdeConfig.silent_) {
		std::cout << "# LearnerSGDE: grid points " << grid->getSize()
				<< std::endl;
	}

	std::shared_ptr<OperationMatrix> C = getRegularizationMatrix(
			configuration.regularizationConfig.regType_);

	for (size_t ref = 0; ref <= configuration.adaptivityConfig.numRefinements_;
			ref++) {
		datadriven::DensitySystemMatrix SMatrix(*grid, samples, *C, 1e-10);
		SMatrix.generateb(rhs);

		if (!configuration.sgdeConfig.silent_) {
			std::cout << "# LearnerSGDE: Solving " << std::endl;
		}

		solver::ConjugateGradients myCG(
				configuration.solverConfig.maxIterations_,
				configuration.solverConfig.eps_);
		myCG.solve(SMatrix, *alpha, rhs, false, false,
				configuration.solverConfig.threshold_);

		if (ref < configuration.adaptivityConfig.numRefinements_) {
			if (!configuration.sgdeConfig.silent_) {
				std::cout << "# LearnerSGDE: Refine grid ... ";
			}

			//Weight surplus with function evaluation at grid points
			OperationEval* opEval = op_factory::createOperationEval(*grid);
			GridIndex* gp;
			DataVector p(numDims);
			DataVector alphaWeight(alpha->getSize());

			for (size_t i = 0; i < gridStorage->size(); i++) {
				gp = gridStorage->get(i);
				gp->getCoords(p);
				alphaWeight[i] = alpha->get(i) * opEval->eval(*alpha, p);
			}

			delete opEval;
			opEval = NULL;

			base::SurplusRefinementFunctor srf(&alphaWeight,
					configuration.adaptivityConfig.noPoints_,
					configuration.adaptivityConfig.threshold_);
			gridGen->refine(&srf);

			if (!configuration.sgdeConfig.silent_) {
				std::cout << "# LearnerSGDE: ref " << ref << "/"
						<< configuration.adaptivityConfig.numRefinements_ - 1
						<< ": " << grid->getStorage()->size() << std::endl;
			}

			alpha->resize(grid->getStorage()->size());
			rhs.resize(grid->getStorage()->size());
			alpha->setAll(0.0);
			rhs.setAll(0.0);
		}
	}
}

} /* namespace datadriven */
} /* namespace SGPP */
