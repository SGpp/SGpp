// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//TODO: Delete this Definition after Refactoring
#define USE_GSL

#ifdef USE_GSL

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/algorithm/ConvergenceMonitor.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>


#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <mpi.h>
#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/application/work_in_progress/MPITaskScheduler.hpp>

#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif


namespace sgpp {
    namespace datadriven {

        using sgpp::base::DataMatrix;
        using sgpp::base::DataVector;

        typedef ClassDensityConntainer ClassDensityContainer;

        /**
 * LearnerSGDEOnOffParallel learns the data using sparse grid density estimation. The
 * system matrix is precomputed and factorized using Eigen-, LU- or
 * Cholesky decomposition (offline step). Then, for each class a density
 * function
 * is computed by solving the system in every iteration (online step).
 * If Cholesky decomposition is chosen, refinement/coarsening can be applied.
 */

        struct LevelIndexPair {
            unsigned long level;
            unsigned long index;
        };

        typedef std::vector<LevelIndexPair> LevelIndexVector;

        struct RefinementResult {
            std::list<LevelIndexVector> addedGridPoints;
            std::list<size_t> deletedGridPointsIndexes;
        };


        class LearnerSGDEOnOffParallel : public LearnerSGDEOnOff {
        public:

            LearnerSGDEOnOffParallel(DBMatDensityConfiguration &dconf, Dataset &trainData, Dataset &testData,
                                     Dataset *validationData, DataVector &classLabels, size_t numClassesInit,
                                     bool usePrior,
                                     double beta, double lambda, MPITaskScheduler &mpiTaskScheduler);

            /**
                                     * Trains the learner with the given dataset.
                                     *
                                     * @param batchSize Size of subset of data points used for each training step
                                     * @param maxDataPasses The number of passes over the whole training data
                                     * @param refinementFunctorType The refinement indicator (surplus, zero-crossings or
                                     * data-based)
                                     * @param refMonitor The refinement strategy (periodic or convergence-based)
                                     * @param refPeriod The refinement interval (if periodic refinement is chosen)
                                     * @param accDeclineThreshold The convergence threshold
                                     *        (if convergence-based refinement is chosen)
                                     * @param accDeclineBufferSize The number of accuracy measurements which are
                                     * used to check
                                     *        convergence (if convergence-based refinement is chosen)
                                     * @param minRefInterval The minimum number of data points (or data batches)
                                     * which have to be
                                     *        processed before next refinement can be scheduled (if
                                     * convergence-based refinement
                                     *        is chosen)
                                     * @param enableCv Specifies whether to perform cross-validation during
                                     * training process or not
                                     * @param nextCvStep Determines when next cross-validation has to be triggered
                                     */
            void train(size_t batchSize, size_t maxDataPasses, std::string refinementFunctorType,
                       std::string refMonitor, size_t refPeriod, double accDeclineThreshold,
                       size_t accDeclineBufferSize, size_t minRefInterval, bool enableCv,
                       size_t nextCvStep);

            /**
             * Trains the learner with the given data batch
             *
             * @param trainData The next data batch to process
             * @param classes The class labels corresponding to the data batch
             * @param doCrossValidation Enable cross-validation
             * @param vectorRefinementResults Vector of pairs containing a list representing indices
             *        of removed grid points and an unsigned int representing added grid
             * points
             */
            void train(Dataset &dataBatch, bool doCrossValidation,
                       std::vector<RefinementResult> *vectorRefinementResults);

            /**
             * Trains the learner with the given data batch that is already split up wrt
             * its different
             * classes.
             *
             * @param trainDataClasses A vector of pairs; Each pair contains the data
             * points that belong to
             *        one class and the corresponding class label
             * @param doCrossValidation Enable cross-validation
             * @param vectorRefinementResults Vector of pairs containing a list representing indices
             * of
             *        removed grid points and an unsigned int representing added grid
             * points
             */
            void train(
                    std::vector<std::pair<base::DataMatrix *, double> > &trainDataClasses,
                    bool doCrossValidation = false,
                    std::vector<RefinementResult> *vectorRefinementResults =
                    nullptr);

            // The final classification error
            double error;


            void updateClassVariablesAfterRefinement(size_t classIndex, RefinementResult *refinementResult,
                                                     DBMatOnlineDE *densEst);

            size_t getDimensionality();

            virtual ~LearnerSGDEOnOffParallel();

            void shutdown();

            void assembleNextBatchData(Dataset *dataBatch, size_t *batchOffset) const;

            void workBatch(Dataset dataset, size_t batchOffset, bool doCrossValidation);

            void mergeAlphaValues(size_t classIndex, size_t gridVersion, DataVector &dataVector, size_t batchSize);

            size_t getCurrentGridVersion(size_t classIndex);

            void setLocalGridVersion(size_t classIndex, size_t gridVersion);

            RefinementResult &getRefinementResult(size_t classIndex);

            void computeNewCholeskyDecomposition(size_t classIndex, size_t gridversion);

            bool checkGridStateConsistent(size_t classIndex);

            static bool isVersionConsistent(size_t version);

        protected:

            std::vector<RefinementResult> *vectorRefinementResults;
            std::vector<size_t> localGridVersions;
            bool workerActive;
            MPITaskScheduler &mpiTaskScheduler;

            size_t
            handleDataAndZeroBasedRefinement(bool preCompute, MultiGridRefinementFunctor *func, size_t idx,
                                             base::Grid &grid,
                                             base::GridGenerator &gridGen) const;


            void allocateClassMatrices(size_t dim, std::vector<std::pair<base::DataMatrix *, double>> &trainDataClasses,
                                       std::map<double, int> &classIndices) const;

            void doRefinementForClass(const std::string &refType, RefinementResult *refinementResult,
                                      const ClassDensityContainer &onlineObjects, bool preCompute,
                                      MultiGridRefinementFunctor *refinementFunctor, size_t classIndex);

            void doRefinementForAll(const std::string &refinementFunctorType,
                                    const std::string &refinementMonitorType,
                                    std::vector<RefinementResult> *vectorRefinementResults,
                                    const ClassDensityContainer &onlineObjects,
                                    ConvergenceMonitor &monitor);


            void printGridSizeStatistics(const char *messageString, ClassDensityContainer &onlineObjects);

            bool checkRefinementNecessary(const std::string &refMonitor, size_t refPeriod, size_t totalInstances,
                                          double currentValidError, double currentTrainError,
                                          size_t numberOfCompletedRefinements,
                                          ConvergenceMonitor &monitor);

            size_t handleSurplusBasedRefinement(DBMatOnlineDE *densEst, Grid &grid, base::GridGenerator &gridGen) const;

            void splitBatchIntoClasses(const Dataset &dataset, size_t dim,
                                       const std::vector<std::pair<DataMatrix *, double>> &trainDataClasses,
                                       std::map<double, int> &classIndices) const;

            size_t assignBatchToWorker(size_t batchOffset, bool doCrossValidation);

            bool checkReadyForRefinement() const;

        };
    }   //namespace datadriven
}  // namespace sgpp

#endif /* USE_GSL */
