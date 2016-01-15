/*
 * DensityRegressionMetaLearner.hpp
 *
 *  Created on: Jan 7, 2016
 *      Author: pfandedd
 */

#include <vector>

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/pde/application/RegularizationConfiguration.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/exception/application_exception.hpp>

namespace SGPP {
namespace datadriven {

class DensityRegressionMetaLearner {
private:
    bool verbose;
    base::DataMatrix &dataset;
    base::DataVector &datasetValues;
    size_t dim;
    base::RegularGridConfiguration gridConfig;
    base::AdpativityConfiguration adaptConfig;
    solver::SLESolverConfiguration solverConfig;
    pde::RegularizationConfiguration regularizationConfig;

    /**
     * generates a regular grid
     * @param dim number of dimensions
     * @return grid
     */
    base::Grid *createRegularGrid(size_t dim);

    float_t calculateMSE(base::Grid &grid, base::DataVector &alpha, base::DataMatrix &testSubset,
            base::DataVector &valuesTestSubset);

    void splitset(std::vector<base::DataMatrix *>& strain, std::vector<base::DataVector *>& strainValues,
            std::vector<base::DataMatrix *>& stest, std::vector<base::DataVector *>& stestValues, size_t kFold,
            bool shuffleDataset, uint32_t shuffleSeed);

    void optimizeLambdaCVGreedy_(size_t kFold, size_t maxLevel, std::vector<base::DataMatrix *> &trainingSets,
            std::vector<base::DataVector *> &trainingSetsValues, std::vector<base::DataMatrix *> &testSets,
            std::vector<base::DataVector *> &testSetsValues, size_t curLevel, float_t lambdaStepSize,
            bool lambdaStepSizeChanged, float_t &bestLambda, float_t &bestMSE);

public:

    DensityRegressionMetaLearner(bool verbose, base::DataMatrix &trainingDataSet, base::DataVector &valuesDataSet,
            base::RegularGridConfiguration gridConfig, base::AdpativityConfiguration adaptConfig,
            solver::SLESolverConfiguration solverConfig, pde::RegularizationConfiguration regularizationConfig);

    float_t optimizeLambdaCV(size_t kFold, float_t lambdaStart, float_t lambdaEnd, size_t lambdaSteps);

    float_t optimizeLambdaCVGreedy(size_t kFold, size_t maxLevel, float_t lambdaStart, float_t lambdaStepSize);

    /**
     * Does the learning step on a given grid, training set and regularization parameter lambda
     *
     * @param grid grid
     * @param train sample set
     * @param trainValues training values
     * @param lambdaReg regularization parameter
     */
    SGPP::base::DataVector train(base::Grid& grid,
    SGPP::base::DataMatrix& train, base::DataVector &trainValues, float_t lambdaReg);
};

}
}
