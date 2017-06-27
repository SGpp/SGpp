// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//TODO: Delete this Definition after Refactoring
#define USE_GSL=1

#ifdef USE_GSL

#ifndef LearnerSGDEOnOffParallel_HPP
#define LearnerSGDEOnOffParallel_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>

#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <mpi.h>

namespace sgpp {
    namespace datadriven {

/**
 * LearnerSGDEOnOffParallel learns the data using sparse grid density estimation. The
 * system matrix is precomputed and factorized using Eigen-, LU- or
 * Cholesky decomposition (offline step). Then, for each class a density
 * function
 * is computed by solving the system in every iteration (online step).
 * If Cholesky decomposition is chosen, refinement/coarsening can be applied.
 */

        struct RefinementResult;
        struct DataBatch;
        struct PendingMPIRequest;

        class LearnerSGDEOnOffParallel : public DBMatOnline {
        public:
            /**
             * Constructor.
             *
             * @param dconf The configuration of the offline object
             * @param trainData The training data
             * @param trainDataLabels The corresponding training labels
             * @param testData The test data
             * @param testDataLabels The corresponding test labels
             * @param validData The validation data
             * @param validDataLabels The corresponding validation labels
             * @param classLabels The class labels (e.g. -1, 1)
             * @param classNumber Total number of classes
             * @param usePrior Determines if prior probabilities should be used to compute
             * class labels
             * @param beta The initial weighting factor
             * @param lambda The initial regularization parameter
             */
            LearnerSGDEOnOffParallel(sgpp::datadriven::DBMatDensityConfiguration &dconf,
                                     base::DataMatrix &trainData,
                                     base::DataVector &trainDataLabels,
                                     base::DataMatrix &testData, base::DataVector &testDataLabels,
                                     base::DataMatrix *validData,
                                     base::DataVector *validDataLabels,
                                     base::DataVector &classLabels, size_t classNumber,
                                     bool usePrior, double beta, double lambda);

            /**
             * Destructor.
             */
            ~LearnerSGDEOnOffParallel();

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
            void train(
                    DataBatch *dataBatch,
                    bool doCrossValidation = false,
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

            /**
             * Returns the accuracy of the classifier measured on the test data.
             *
             * @return The classification accuracy measured on the test data
             */
            double getAccuracy();

            /**
             * Predicts the class labels of the test data points.
             *
             * @param The test data points
             * @return A vector containing the predicted class labels
             */
            sgpp::base::DataVector predict(base::DataMatrix &test);

            /**
             * Predicts the class label of the given data point.
             *
             * @param p The data point
             * @return The predicted class label
             */
            int predict(base::DataVector &p);

            /**
             * Error evaluation required for convergence-based refinement.
             *
             * @param data The data points to measure the error on
             * @param labels The corresponding class labels
             * @param errorType The error type (only "Acc" possible, i.e. classification
             * error
             *        based on accuracy)
             * @return The error evaluation
             */
            double getError(base::DataMatrix &data, base::DataVector &labels,
                            std::string errorType);

            /**
             * Stores classified data, grids and density function evaluations to csv
             * files.
             */
            void storeResults();

            /**
             * Returns the values of all density functions for a specified data point.
             *
             * @param point The point for which the density functions should be evaluated
             * @return The function evaluations
             */
            base::DataVector getDensities(base::DataVector &point);

            /**
             * Sets the cross-validation parameters.
             * They get directly passed to the DBMatOnlineDE class-instance.
             *
             * @param lambdaStep Defines how many different lambdas are tried out
             * @param lambdaStart The smallest possible lambda
             * @param lambdaEnd The biggest possible lambda
             * @param test The test matrix
             * @param testRes The results of the points in the test matrix
             * @param logscale Indicates whether the values between lambdaStart
             *        and lambdaEnd are searched using logscale or not
             */
            void setCrossValidationParameters(int lambdaStep, double lambdaStart,
                                              double lambdaEnd, base::DataMatrix *test,
                                              base::DataMatrix *testRes, bool logscale);

            /**
            * In case of crossvalidation, returns the current best lambda.
            *
            * @return The lambda value
            */
            // double getBestLambda();

            /**
            * Initialization of online objects in case of Eigen- or LU-decomposition.
            */
            void init();

            /**
            * Returns the number of existing classes.
            *
            * @return The number of classes
            */
            size_t getNumClasses();

            /**
             * Returns the density functions mapped to class labels.
             *
             * @return The density function objects mapped to class labels
             */
            std::vector<std::pair<DBMatOnlineDE *, double> > *getDestFunctions();

            // Stores prior values mapped to class labels
            std::map<double, double> prior;

            // The final classification error
            double error;

            // A vector to store error evaluations
            base::DataVector avgErrors;

        protected:
            // The training data
            base::DataMatrix &trainData;
            // The corresponding training class labels
            base::DataVector &trainLabels;
            // The test data
            base::DataMatrix &testData;
            // The corresponding test class labels
            base::DataVector &testLabels;
            // The validation data
            base::DataMatrix *validData;
            // The corresponding validation class labels
            base::DataVector *validLabels;

            // The class labels (e.g -1, 1)
            base::DataVector classLabels;
            // The total number of different classes
            size_t classNumber;

            // Indicates whether the model has been trained or not
            bool trained;
            // Indicates whether the learner has been initialized or not
            bool initDone;
            // Specifies whether prior should be used for class prediction or not
            bool usePrior;
            // Weighting factor
            double beta;

            // The offline object (contains decomposed matrix)
            DBMatOffline *offline;
            // The online objects (density functions)
            std::vector<std::pair<DBMatOnlineDE *, double> > *destFunctions;

            // Counter for total number of data points processed within ona data pass
            size_t processedPoints;

            // Cross-validation parameters
            int cvSaveLambdaStep;
            double cvSaveLambdaStart;
            double cvSaveLambdaEnd;
            bool cvSaveLogscale;
            bool cvSaved;
            base::DataMatrix *cvSaveTest;
            base::DataMatrix *cvSaveTestRes;

            //Pending MPI Requests
            std::vector<PendingMPIRequest> pendingMPIRequests;
            int mpiWorldSize;

            //Batches assigned by master
            std::vector<DataBatch> assignedBatches;

            size_t
            handleSurplusBasedRefinement(DBMatOnlineDE *densEst, base::Grid *grid,
                                         base::GridGenerator &gridGen) const;

            size_t
            handleDataAndZeroBasedRefinement(bool preCompute, MultiGridRefinementFunctor *func, size_t idx,
                                             base::Grid *grid,
                                             base::GridGenerator &gridGen) const;

            bool checkRefinementNecessary(const std::string &refMonitor, size_t refPeriod,
                                          size_t totalInstances, double currentValidError,
                                          double currentTrainError, size_t numberOfCompletedRefinements,
                                          std::shared_ptr<ConvergenceMonitor> &monitor);

            void assembleNextBatchData(size_t batchSize, DataBatch *dataBatch,
                                       size_t dataDimensionality,
                                       size_t *batchOffset) const;

            void
            printGridSizeStatistics(std::vector<std::pair<DBMatOnlineDE *, double>> *onlineObjects,
                                    const char *messageString);

            void allocateClassMatrices(size_t dim, std::vector<std::pair<base::DataMatrix *, double>> &trainDataClasses,
                                       std::map<double, int> &classIndices) const;

            void
            splitBatchIntoClasses(const DataBatch *dataBatch, size_t dim,
                                  const std::vector<std::pair<base::DataMatrix *, double>> &trainDataClasses,
                                  std::map<double, int> &classIndices) const;

            void initMPI();

            void synchronizeEndOfDataPass();

            bool isMaster() const;

            void sendGridComponentsUpdate(std::vector<RefinementResult> *refinementResults);


            void waitForMPIRequestsToComplete();

            void processCompletedMPIRequests();

            void sendCommandIDToWorkers(sgpp::datadriven::MPI_COMMAND_ID commandID) const;

            void processIncomingMPICommands(sgpp::datadriven::MPI_Packet *mpiPacket);

            void receiveGridComponentsUpdate(sgpp::datadriven::RefinementResultNetworkMessage *networkMessage);

            void updateVariablesAfterRefinement(RefinementResult *refinementResult, size_t classIndex,
                                                DBMatOnlineDE *densEst) const;

            void doRefinementForClass(const std::string &refType, RefinementResult *refinementResult,
                                      const std::pair<DBMatOnlineDE *, double> *onlineObjects, bool preCompute,
                                      MultiGridRefinementFunctor *refinementFunctor, size_t classIndex) const;

            void doRefinementForAll(const std::string &refinementFunctorType,
                                    const std::string &refinementMonitorType,
                                    std::vector<RefinementResult> *vectorRefinementResults,
                                    const std::vector<std::pair<DBMatOnlineDE *, double>> *onlineObjects,
                                    std::shared_ptr<ConvergenceMonitor> &monitor);
        };
    }   //namespace datadriven
}  // namespace sgpp

#endif /* LearnerSGDEOnOffParallel_HPP */

#endif /* USE_GSL */
