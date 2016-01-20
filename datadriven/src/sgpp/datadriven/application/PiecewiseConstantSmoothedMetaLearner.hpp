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

class PiecewiseConstantSmoothedRegressionMetaLearner {
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

//    void splitset(std::vector<base::DataMatrix *>& strain, std::vector<base::DataVector *>& strainValues,
//            std::vector<base::DataMatrix *>& stest, std::vector<base::DataVector *>& stestValues, size_t kFold,
//            bool shuffleDataset, uint32_t shuffleSeed);

    void optimizeLambdaLog_(size_t kFold, size_t maxLevel, float_t fastApproximationMSE,
            size_t fastApproximationMaxLevel, std::vector<base::DataMatrix> &trainingSets,
            std::vector<base::DataVector> &trainingSetsValues, std::vector<base::DataMatrix> &testSets,
            std::vector<base::DataVector> &testSetsValues, size_t curLevel, float_t lambdaLogStepSize,
            float_t &bestLogLambda, float_t &bestMSE);

public:

    PiecewiseConstantSmoothedRegressionMetaLearner(bool verbose, base::DataMatrix &trainingDataSet, base::DataVector &valuesDataSet,
            base::RegularGridConfiguration gridConfig, base::AdpativityConfiguration adaptConfig,
            solver::SLESolverConfiguration solverConfig, pde::RegularizationConfiguration regularizationConfig);

    float_t optimizeLambdaLog(size_t kFold, size_t maxLevel, float_t fastApproximationMSE,
            size_t fastApproximationMaxLevel);

    void optimizeLambdaLog(size_t kFold, size_t maxLevel, float_t fastApproximationMSE,
            size_t fastApproximationMaxLevel, std::shared_ptr<base::Grid> &bestGrid,
            std::shared_ptr<base::DataVector> &bestAlpha, float_t &lambdaOpt);

    /**
     * Does the learning step on a given grid, training set and regularization parameter lambda
     *
     * @param train sample set
     * @param trainValues training values
     * @param lambdaReg regularization parameter
     * @param fastApproximationMSE mse for stopping piecewise constant approximation tree creation
     * @param fastApproximationMSE maximum level for stopping piecewise constant approximation tree creation
     * @param grid grid
     * @param grid grid coefficients
     */
    void train(base::DataMatrix &train, base::DataVector &trainValues, float_t lambda, float_t fastApproximationMSE,
            size_t fastApproximationMaxLevel, std::shared_ptr<base::Grid> &grid,
            std::shared_ptr<base::DataVector> &alpha);

    float_t calculateMSE(base::Grid &grid, base::DataVector &alpha, base::DataMatrix &testSubset,
            base::DataVector &valuesTestSubset, bool verbose = false);
};

}
}
