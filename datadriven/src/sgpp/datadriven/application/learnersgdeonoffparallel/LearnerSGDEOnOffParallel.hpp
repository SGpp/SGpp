// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <mpi.h>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPITaskScheduler.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/AuxiliaryStructures.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/RefinementHandler.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorConvergence.hpp>

#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

/**
* LearnerSGDEOnOffParallel learns the data using sparse grid density estimation.
* The
* system matrix is precomputed and factorized using Eigen-, LU- or
* Cholesky decomposition (offline step). Then, for each class a density
* function is computed by solving the system in every iteration (online step).
* If Cholesky decomposition is chosen, refinement/coarsening can be applied.
* This learner uses MPI to parallelize the learning phase across multiple nodes.
*/

class LearnerSGDEOnOffParallel {
 public:
  LearnerSGDEOnOffParallel(
      sgpp::base::RegularGridConfiguration &gridConfig,
      sgpp::base::AdaptivityConfiguration &adaptivityConfig,
      sgpp::datadriven::RegularizationConfiguration &regularizationConfig,
      sgpp::datadriven::DensityEstimationConfiguration &densityEstimationConfig, Dataset &trainData,
      Dataset &testData, Dataset *validationData, DataVector &classLabels, size_t numClassesInit,
      bool usePrior, double beta, MPITaskScheduler &mpiTaskScheduler);

  /**
   * Trains the learner with the given dataset.
   *
   * @param batchSize Size of subset of data points used for each training step
   * @param maxDataPasses The number of passes over the whole training data
   * @param refinementFunctorType The refinement indicator (surplus,
   * zero-crossings or
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
   */
  void trainParallel(size_t batchSize, size_t maxDataPasses, std::string refinementFunctorType,
                     std::string refMonitor, size_t refPeriod, double accDeclineThreshold,
                     size_t accDeclineBufferSize, size_t minRefInterval);

  /**
   * Trains the learner with the given data batch
   *
   * @param dataBatch The next data batch to process
   * @param doCrossValidation Enable cross-validation
   */
  void train(Dataset &dataBatch, bool doCrossValidation);

  /**
   * Trains the learner with the given data batch that is already split up wrt
   * its different classes.
   *
   * @param trainDataClasses A vector of pairs; Each pair contains the data
   * points that belong to one class and the corresponding class label
   * @param doCrossValidation Enable cross-validation
   *
   */
  void train(std::vector<std::pair<DataMatrix *, double>> &trainDataClasses,
             bool doCrossValidation);

  /**
   * Returns the dimensionality of the learner as determined from its training
   * set
   *
   * @return The data dimensionality
   */
  size_t getDimensionality();

  /**
   * Runs MPI finalize when destructing the learner
   */
  virtual ~LearnerSGDEOnOffParallel();

  /**
   * If this is run on master, it issues shutdown requests to all workers and
   * waits
   * for them to return.
   * If this is run on a worker, it sets the shutdown flag.
   */
  void shutdownMPINodes();

  /**
   * Copies the data from the training set into the data batch
   *
   * @param dataBatch Batch of data to fill, with set dimensionality and size
   * @param batchOffset The offset in the training data from which to start
   * copying
   */
  void assembleNextBatchData(Dataset *dataBatch, size_t *batchOffset) const;

  /**
   * Train from a batch. Will wait until all grids are consistent, fill the
   * dataset,
   * learn from the dataset and send the new alpha vector to the master
   *
   * @param dataset An empty dataset with size and dimension set.
   * @param batchOffset The offset from the start of the training set to
   * assemble the batch from.
   * @param doCrossValidation Whether to cross validate results.
   */
  void workBatch(Dataset dataset, size_t batchOffset, bool doCrossValidation);

  /**
   * Merge alpha values received from a remote process into the local alpha
   * vector.
   *
   * @param classIndex The class to which the alpha vector belongs
   * @param remoteGridVersion The remote grid version this alpha vector was
   * trained on
   * @param dataVector The alpha vector itself
   * @param batchOffset The offset from the start of the training set this
   * vector was trained from
   * @param batchSize The size of the batch this vector was trained from
   * @param isLastPacketInSeries Whether this merge is the last merge in several
   * for the same class and batch
   */
  void mergeAlphaValues(size_t classIndex, size_t remoteGridVersion, DataVector dataVector,
                        size_t batchOffset, size_t batchSize, bool isLastPacketInSeries);

  /**
   * Returns the internally stored current version of the grid
   *
   * @param classIndex The class of the grid to search for
   * @return The current version of the grid
   */
  size_t getLocalGridVersion(size_t classIndex);

  /**
   * Set the grid version
   *
   * @param classIndex The class of the grid to search for
   * @param gridVersion The new version of the grid
   */
  void setLocalGridVersion(size_t classIndex, size_t gridVersion);

  /**
   * Update the system matrix decomposition after a refinement step.
   * This will wait for the receiving of refinement results to complete.
   * After computation, the system matrix is sent back to the master
   *
   * @param classIndex The class for which to update the system matrix
   * decomposition
   * @param gridVersion The new grid version to set after updating the matrix
   */
  void computeNewSystemMatrixDecomposition(size_t classIndex, size_t gridVersion);

  /**
   * Check whether the grid is in a final state where learning can occur.
   * This is not the case while receiving refinement results or updating
   * the system matrix decomposition.
   *
   * @param classIndex The class for which to check consistency.
   * @return Whether the grid is currently in a consistent state
   */
  bool checkGridStateConsistent(size_t classIndex);

  /**
   * Check whether a specific grid version is consistent, i.e. whether it is
   * higher than MINIMUM_CONSISTENT_GRID_VERSION
   *
   * @param version The version of the grid to check against
   * @return Whether the version indicates consistency.
   */
  static bool isVersionConsistent(size_t version);

  /**
   * Gets a reference to the currently installed MPI Scheduler.
   * The scheduler assigns tasks of variable or static size to workers.
   *
   * @return A reference to the installed MPI Task Scheduler
   */
  MPITaskScheduler &getScheduler();

  /**
   * Gets the DBMatOffline object
   *
   * @return The DBMatOffline object
   */
  std::unique_ptr<DBMatOffline> &getOffline();

  /**
   * Check whether all grids are not in a temporarily inconsistent state.
   *
   * @return Whether all grids are consistent
   */
  bool checkAllGridsConsistent();

  /**
   * Returns a reference to the currently used training data set
   *
   * @return A reference to the training data set
   */
  Dataset &getTrainData();

  /**
   * Returns a reference to the currently used test data set
   *
   * @return A reference to the test data set
   */
  Dataset *getValidationData();

  /**
   * Returns a reference to the refinement handler, that contains logic to
   * handle
   * the master's refinement cycles
   *
   * @return A reference to the refinement handler
   */
  RefinementHandler &getRefinementHandler();

  /**
   * Asks the scheduler where to assign the next batch to and sends the MPI
   * request.
   *
   * @param batchOffset Starting offset of the new batch
   * @param doCrossValidation Whether the client should do cross-validation
   * @return The size of the batch assigned by the scheduler
   */
  size_t assignBatchToWorker(size_t batchOffset, bool doCrossValidation);

  /**
   * Retrieves the grid for a certain class
   * @param classIndex the index of the desired class
   * @return the underlying grid
   */
  Grid &getGrid(size_t classIndex);

  /**
  * Returns the number of existing classes.
  *
  * @return The number of classes
  */
  size_t getNumClasses() const;

  /**
   * Returns the accuracy of the classifier measured on the test data.
   *
   * @return The classification accuracy measured on the test data
   */
  double getAccuracy() const;

  /**
   * Predicts the class labels of the test data points.
   *
   * @param test The data points for which labels will be precicted
   * @param classLabels vector containing the predicted class labels
   */
  void predict(DataMatrix &test, DataVector &classLabels) const;

  /**
   * Error evaluation required for convergence-based refinement.
   *
   * @param dataset The data to measure the error on
   * @return The error evaluation
   */
  double getError(Dataset &dataset) const;

  /**
   * Updates the surplus vector of a certain class
   *
   * @param classIndex the index of the class
   * @param deletedPoints a list of indexes of deleted points (coarsening)
   * @param newPoints the number of new grid points (refinemenet)
   */
  void updateAlpha(size_t classIndex, std::list<size_t> *deletedPoints, size_t newPoints);

  /**
   * Returns the density functions mapped to class labels.
   *
   * @return The density function objects mapped to class labels
   */
  std::vector<std::pair<std::unique_ptr<DBMatOnlineDE>, size_t>> &getDensityFunctions();

 protected:
  // Grids TODO(fuchsgruber): Move outwards (just in this class so that it
  // compiles...)
  std::vector<std::unique_ptr<Grid>> grids;
  // Surplusses TODO(fuchsgruber): Move alphas outwards (just in this class so
  // that it compiles)
  std::vector<DataVector *> alphas;

  // The training data
  Dataset &trainData;
  // The test data
  Dataset &testData;

  // The (optional) validationData
  Dataset *validationData;

  // The class labels (e.g -1, 1)
  DataVector classLabels;
  // The total number of different classes
  size_t numClasses;
  // Specifies whether prior should be used for class prediction or not
  bool usePrior;
  // Stores prior values mapped to class labels
  std::map<double, double> prior;
  // Weighting factor
  double beta;
  // Indicates whether the model has been trained or not
  bool trained;

  // Contains the offline object that was cloned into all other classes
  std::unique_ptr<DBMatOffline> offline;
  // Contains all offline objects
  std::vector<std::unique_ptr<DBMatOffline>> offlineContainer;
  // The online objects (density functions)
  std::vector<std::pair<std::unique_ptr<DBMatOnlineDE>, size_t>> densityFunctions;

  // Counter for total number of data points processed within ona data pass
  size_t processedPoints;

  // The final classification error
  // double error;

  // A vector to store error evaluations
  DataVector avgErrors;

  /**
   * Vector that holds the grid version for every class
   */
  std::vector<size_t> localGridVersions;

  /**
   * Boolean used to detect when a shutdown of a worker has been requested
   */
  bool workerActive;

  /**
   * Reference to the currently installed MPI Task Scheduler
   */
  MPITaskScheduler &mpiTaskScheduler;

  /**
   * Instance of the currently installed refinement handler
   */
  RefinementHandler refinementHandler;

  // Configuration for the grid TODO(fuchsgruber): Move outwards
  sgpp::base::GeneralGridConfiguration &gridConfig;

  // Configuration for the adaptivity TODO(fuchsgruber): Move outwards
  sgpp::base::AdaptivityConfiguration &adaptivityConfig;

  // Configuration for the regularization TODO(fuchsgruber): Move outwards
  sgpp::datadriven::RegularizationConfiguration &regularizationConfig;

  // Configuration for the density estimation TODO(fuchsgruber): Move outwards
  sgpp::datadriven::DensityEstimationConfiguration &densityEstimationConfig;

  /**
   * Allocates memory for every class to hold training data before learning
   *
   * @param dim The dimensionality of the current problem
   * @param trainDataClasses Storage that will be allocated that holds space for
   * data and label
   * @param classIndices A map of each classes label to its index
   */
  void allocateClassMatrices(size_t dim,
                             std::vector<std::pair<DataMatrix *, double>> &trainDataClasses,
                             std::map<double, int> &classIndices) const;

  /**
   * Do an entire refinement cycle for all classes.
   *
   * @param refinementFunctorType String constant specifying the functor to use
   * in refinement
   * @param refinementMonitorType String constant specifying the monitor to use
   * in refinement
   * @param onlineObjects Reference to the online objects for density estimation
   * @param monitor The setup of the convergence monitor for refinement
   */
  void doRefinementForAll(
      const std::string &refinementFunctorType, const std::string &refinementMonitorType,
      const std::vector<std::pair<std::unique_ptr<DBMatOnlineDE>, size_t>> &onlineObjects,
      RefinementMonitor &monitor);

  /**
   * Shows grid size statistics along with a message
   *
   * @param messageString The message to display alongside the statistics
   * @param onlineObjects The current density estimation objects
   */
  void printGridSizeStatistics(
      const char *messageString,
      std::vector<std::pair<std::unique_ptr<DBMatOnlineDE>, size_t>> &onlineObjects);

  void splitBatchIntoClasses(const Dataset &dataset, size_t dim,
                             const std::vector<std::pair<DataMatrix *, double>> &trainDataClasses,
                             std::map<double, int> &classIndices) const;

  /**
   * Wait for all grids to reach a consistent state before continuing
   */
  void waitForAllGridsConsistent();
};
}  // namespace datadriven
}  // namespace sgpp
