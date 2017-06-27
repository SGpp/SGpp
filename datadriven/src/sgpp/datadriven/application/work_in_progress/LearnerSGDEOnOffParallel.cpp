// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//TODO: Delete this Definition after Refactoring
#define USE_GSL=1


#ifdef USE_GSL

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/algorithm/ConvergenceMonitor.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/application/work_in_progress/LearnerSGDEOnOffParallel.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>

#include <sgpp/parallel/tools/MPI/GlobalMPIComm.hpp>
#include <sgpp/parallel/tools/MPI/MPICommunicator.hpp>

#include <ctime>

namespace sgpp {
    namespace datadriven {

        static const int MPI_PACKET_MAX_PAYLOAD_SIZE = 512;

        struct RefinementResult {
            size_t numAddedGridPoints;
            std::list<size_t> deletedGridPointsIndexes;
        };

        struct DataBatch {
            // pointer to the next batch (data points + class labels) to be processed
            sgpp::base::DataMatrix *dataPoints;
            sgpp::base::DataVector *classLabels;
        };

        struct PendingMPIRequest {
            MPI_Request request;
            MPI_Packet *buffer;
            std::function<void(void *)> callback;
            bool disposeAfterCallback;
        };

        struct MPI_Packet {
            MPI_COMMAND_ID commandID;
            unsigned char payload[MPI_PACKET_MAX_PAYLOAD_SIZE];
        };

        enum MPI_COMMAND_ID {
            UPDATE_GRID,
            MERGE_GRID,
            ASSIGN_BATCH
        };

        enum MPI_COMMAND_TAG {

        };

        const int MPI_MASTER_RANK = 0;


        LearnerSGDEOnOffParallel::LearnerSGDEOnOffParallel(
                sgpp::datadriven::DBMatDensityConfiguration &dconf,
                sgpp::base::DataMatrix &trainData, sgpp::base::DataVector &trainDataLabels,
                sgpp::base::DataMatrix &testData, sgpp::base::DataVector &testDataLabels,
                sgpp::base::DataMatrix *validData, sgpp::base::DataVector *validDataLabels,
                sgpp::base::DataVector &classLabels, size_t classNumber,
                bool usePrior, double beta, double lambda)
                : trainData(trainData),
                  trainLabels(trainDataLabels),
                  testData(testData),
                  testLabels(testDataLabels),
                  validData(validData),
                  validLabels(validDataLabels),
                  classLabels(classLabels),
                  classNumber(classNumber),
                  trained(false),
                  initDone(false),
                  usePrior(usePrior),
                  beta(beta),
                  destFunctions(nullptr) {
            offline = new DBMatOffline(dconf);
            offline->buildMatrix();
            // clock_t begin = clock();
            offline->decomposeMatrix();
            // clock_t end = clock();
            // double elapsed_secs = double(end-begin)/CLOCKS_PER_SEC;
            // std::cout << "#Decompose matrix: " << elapsed_secs << std::endl;

            readOffline(offline);  // set offlineObject_ of DBMatOnline class

            init();

            cvSaved = false;
            processedPoints = 0;

            for (size_t i = 0; i < classNumber; i++) {
                prior.insert(std::pair<double, double>(classLabels[i], 0.0));
            }

            setLambda(lambda);
        }

        void LearnerSGDEOnOffParallel::init() {
            if (initDone) {
                return;
            }

            // if the Cholesky decomposition is chosen declare separate Online-objects for
            // every class
            if (offline->getConfig()->decomp_type_ == DBMatDecompChol) {
                destFunctions = new std::vector<std::pair<DBMatOnlineDE *, double> >(classNumber);
                // every class gets his own online object
                for (size_t i = 0; i < classNumber; i++) {
                    DBMatOnlineDE *densEst = new DBMatOnlineDE(beta);
                    //This is the difference
                    DBMatOffline *offlineRead = new DBMatOffline(*offline);
                    densEst->readOffline(offlineRead);
                    std::pair<DBMatOnlineDE *, double> pdest(densEst, classLabels[i]);
                    (*destFunctions)[i] = pdest;
                }
            } else {
                destFunctions = new std::vector<std::pair<DBMatOnlineDE *, double> >(classNumber);
                for (size_t idx = 0; idx < classNumber; idx++) {
                    DBMatOnlineDE *densEst = new DBMatOnlineDE(beta);
                    densEst->readOffline(offline);
                    std::pair<DBMatOnlineDE *, double> pdest(densEst, classLabels[idx]);
                    (*destFunctions)[idx] = pdest;
                }
                if (cvSaved) {
                    setCrossValidationParameters(cvSaveLambdaStep, cvSaveLambdaStart,
                                                 cvSaveLambdaEnd, cvSaveTest, cvSaveTestRes,
                                                 cvSaveLogscale);
                }
            }
            initMPI();
            initDone = true;
        }

        LearnerSGDEOnOffParallel::~LearnerSGDEOnOffParallel() {
            if (destFunctions != nullptr) {
                for (size_t i = 0; i < destFunctions->size(); i++) {
                    delete (*destFunctions)[i].first;
                }
                delete destFunctions;
            }
        }

        void LearnerSGDEOnOffParallel::train(size_t batchSize, size_t maxDataPasses,
                                             std::string refinementFunctorType, std::string refMonitor,
                                             size_t refPeriod, double accDeclineThreshold,
                                             size_t accDeclineBufferSize, size_t minRefInterval,
                                             bool enableCv, size_t nextCvStep) {
            // counts total number of processed data points
            size_t numProcessedDataPoints = 0;
            // pointer to the next batch (data points + class labels) to be processed
            DataBatch dataBatch;

            // contains list of removed grid points and number of added grid points
            // (is updated in each refinement/coarsening step)
            std::vector<RefinementResult> *vectorRefinementResults =
                    new std::vector<RefinementResult>(classNumber);

            // initialize counter for dataset passes
            size_t completedDataPasses = 0;

            // initialize refinement variables
            double currentValidError = 0.0;
            double currentTrainError = 0.0;
            // create convergence monitor object
            std::shared_ptr<ConvergenceMonitor> monitor(new ConvergenceMonitor(
                    accDeclineThreshold, accDeclineBufferSize, minRefInterval));

            // counts number of performed refinements
            size_t numberOfCompletedRefinements = 0;

            // coarsening
            // size_t coarseCnt = 0;
            // size_t maxCoarseNum = 5;
            // size_t coarsePeriod = 50;
            // size_t coarseNumPoints = 1;
            // double coarseThreshold = 1.0;

            std::vector<std::pair<DBMatOnlineDE *, double> > *onlineObjects;

            size_t dim = trainData.getNcols();
            // determine number of batches to process
            size_t numBatch = trainData.getNrows() / batchSize;

            onlineObjects = getDestFunctions();

            // print initial grid size
            printGridSizeStatistics(onlineObjects, "#Initial grid size of grid ");

            // auxiliary variable for accuracy (error) measurement
            double acc = getAccuracy();

            avgErrors.append(1.0 - acc);

            // main loop which performs the training process
            while (completedDataPasses < maxDataPasses) {
                std::cout << "#batch-size: " << batchSize << std::endl;
                std::cout << "#batches to process: " << numBatch << std::endl;

                // data point counter - determines offset when selecting next batch
                size_t cnt = 0;

                // iterate over total number of batches
                for (size_t currentBatchNum = 1; currentBatchNum <= numBatch; currentBatchNum++) {
                    // check if cross-validation should be performed
                    bool doCrossValidation = false;
                    if (enableCv) {
                        if (nextCvStep == currentBatchNum) {
                            doCrossValidation = true;
                            nextCvStep *= 5;
                        }
                    }
                    // assemble next batch
                    assembleNextBatchData(batchSize, &dataBatch, dim, &cnt);

                    // train the model with current batch
                    train(&dataBatch, doCrossValidation, vectorRefinementResults);

                    numProcessedDataPoints += (dataBatch.dataPoints)->getNrows();

                    // access DBMatOnlineDE-objects of all classes in order
                    // to apply adaptivity to the specific sparse grids later on
                    onlineObjects = getDestFunctions();

                    // Refinement only occurs on the Master Node
                    if (isMaster()) {

                        // check if refinement should be performed
                        if (checkRefinementNecessary(refMonitor, refPeriod, numProcessedDataPoints, currentValidError,
                                                     currentTrainError, numberOfCompletedRefinements, monitor)) {
                            // if the Cholesky decomposition is chosen as factorization method
                            // refinement
                            // and coarsening methods can be applied

                            std::cout << "refinement at iteration: " << numProcessedDataPoints << std::endl;
                            doRefinementForAll(refinementFunctorType, refMonitor, vectorRefinementResults,
                                               onlineObjects, monitor);
                            numberOfCompletedRefinements += 1;

                            //TODO If not master, the grid needs to be adjusted here
                            if (isMaster()) {
                                //TODO Send and Receive Delta
                                //TODO Adjust Grid
                                sendGridComponentsUpdate(vectorRefinementResults);
                            }
                        }


                    } else {
                        //TODO: Otherwise, merge grids if any are available
                    }
                    delete dataBatch.dataPoints;
                    delete dataBatch.classLabels;

                    // save current error
                    if (numProcessedDataPoints % 10 == 0) {
                        acc = getAccuracy();
                        avgErrors.append(1.0 - acc);
                    }
                }

                // Synchronize end of Data Pass
                synchronizeEndOfDataPass();

                completedDataPasses++;
                processedPoints = 0;
            }  // end while

            std::cout << "#Training finished" << std::endl;

            error = 1.0 - getAccuracy();

            // delete offline;
            delete vectorRefinementResults;
        }

        bool LearnerSGDEOnOffParallel::isMaster() const {
            return sgpp::parallel::myGlobalMPIComm->getMyRank() == MPI_MASTER_RANK;
        }

        void LearnerSGDEOnOffParallel::printGridSizeStatistics(
                std::vector<std::pair<DBMatOnlineDE *, double>> *onlineObjects,
                const char *messageString) {
            for (size_t i = 0; i < classNumber; i++) {
                DBMatOnlineDE *densEst = (*onlineObjects)[i].first;
                base::Grid *grid = densEst->getOffline()->getGridPointer();
                std::cout << messageString << i << " : " << grid->getSize()
                          << std::endl;
            }
        }

        void LearnerSGDEOnOffParallel::doRefinementForAll(const std::string &refinementFunctorType,
                                                          const std::string &refinementMonitorType,
                                                          std::vector<RefinementResult> *vectorRefinementResults,
                                                          const std::vector<std::pair<DBMatOnlineDE *, double>> *onlineObjects,
                                                          std::shared_ptr<ConvergenceMonitor> &monitor) {
            // acc = getAccuracy();
            // avgErrors.append(1.0 - acc);

            // bundle grids and surplus vector pointer needed for refinement
            // (for zero-crossings refinement, data-based refinement)
            std::vector<sgpp::base::Grid *> grids;
            std::vector<sgpp::base::DataVector *> alphas;
            for (size_t i = 0; i < this->getNumClasses(); i++) {
                DBMatOnlineDE *densEst = (*onlineObjects)[i].first;
                grids.push_back(densEst->getOffline()->getGridPointer());
                alphas.push_back(densEst->getAlpha());
            }
            bool levelPenalize =
                    false;               // Multiplies penalzing term for fine levels
            bool preCompute = true;  // Precomputes and caches evals for zrcr
            sgpp::datadriven::MultiGridRefinementFunctor *func = nullptr;

            // Zero-crossing-based refinement
            sgpp::datadriven::ZeroCrossingRefinementFunctor funcZrcr =
                    *(new sgpp::datadriven::ZeroCrossingRefinementFunctor(
                            grids, alphas, this->offline->getConfig()->ref_noPoints_,
                            levelPenalize, preCompute));

            // Data-based refinement. Needs a problem dependent coeffA. The values
            // can be determined by testing (aim at ~10 % of the training data is
            // to be marked relevant). Cross-validation or similar can/should be
            // employed
            // to determine this value.
            //TODO: Investigate this
            std::vector<double> coeffA;
            coeffA.push_back(1.2);  // ripley 1.2
            coeffA.push_back(1.2);  // ripley 1.2
            base::DataMatrix *trainDataRef = &this->trainData;
            base::DataVector *trainLabelsRef = &this->trainLabels;
            sgpp::datadriven::DataBasedRefinementFunctor funcData =
                    *(new sgpp::datadriven::DataBasedRefinementFunctor(
                            grids, alphas, trainDataRef, trainLabelsRef,
                            this->offline->getConfig()->ref_noPoints_, levelPenalize, coeffA));
            if (refinementFunctorType == "zero") {
                func = &funcZrcr;
            } else if (refinementFunctorType == "data") {
                func = &funcData;
            }

            // perform refinement/coarsening for each grid
            for (size_t idx = 0; idx < this->getNumClasses(); idx++) {
                this->doRefinementForClass(refinementFunctorType,
                                           &((*vectorRefinementResults)[idx]),
                                           &((*onlineObjects)[idx]),
                                           preCompute, func, idx);
            }
            if (refinementMonitorType == "convergence") {
                monitor->nextRefCnt = monitor->minRefInterval;
            }
        }

        void LearnerSGDEOnOffParallel::doRefinementForClass(const std::string &refType,
                                                            RefinementResult *refinementResult,
                                                            const std::pair<DBMatOnlineDE *, double> *onlineObjects,
                                                            bool preCompute,
                                                            MultiGridRefinementFunctor *refinementFunctor,
                                                            size_t classIndex) const {
            // perform refinement/coarsening for grid which corresponds to current
            // index
            std::cout << "Refinement and coarsening for class: " << classIndex
                      << std::endl;
            DBMatOnlineDE *densEst = (*onlineObjects).first;
            base::Grid *grid = densEst->getOffline()->getGridPointer();
            std::cout << "Size before adaptivity: " << grid->getSize()
                      << std::endl;

            base::GridGenerator &gridGen = grid->getGenerator();

            size_t numberOfNewPoints = 0;

            if (refType == "surplus") {
                numberOfNewPoints = handleSurplusBasedRefinement(densEst, grid, gridGen);

            } else if ((refType == "data") || (refType == "zero")) {
                numberOfNewPoints = handleDataAndZeroBasedRefinement(preCompute, refinementFunctor, classIndex, grid,
                                                                     gridGen);
            }

            std::cout << "grid size after adaptivity: " << grid->getSize()
                      << std::endl;

            refinementResult->numAddedGridPoints = numberOfNewPoints;

            updateVariablesAfterRefinement(refinementResult, classIndex, densEst);
        }

        void LearnerSGDEOnOffParallel::updateVariablesAfterRefinement(
                RefinementResult *refinementResult,
                size_t classIndex,
                DBMatOnlineDE *densEst) const {

            // apply grid changes to the Cholesky factorization
            densEst->getOffline()->choleskyModification(
                    refinementResult->numAddedGridPoints,
                    refinementResult->deletedGridPointsIndexes,
                    densEst->getBestLambda());
            // update alpha vector
            densEst->updateAlpha(
                    &refinementResult->deletedGridPointsIndexes,
                    refinementResult->numAddedGridPoints);
        }

        void LearnerSGDEOnOffParallel::assembleNextBatchData(size_t batchSize,
                                                             DataBatch *dataBatch,
                                                             size_t dataDimensionality,
                                                             size_t *batchOffset) const {
            base::DataMatrix *batch =
                    new base::DataMatrix(batchSize, dataDimensionality);
            base::DataVector *batchLabels =
                    new base::DataVector(batchSize);
            for (size_t j = 0; j < batchSize; j++) {
                base::DataVector dataPoint(dataDimensionality);
                trainData.getRow(j + *batchOffset, dataPoint);
                double y = trainLabels.get(j + *batchOffset);
                batch->setRow(j, dataPoint);
                batchLabels->set(j, y);
            }
            dataBatch->dataPoints = batch;
            dataBatch->classLabels = batchLabels;

            *batchOffset += batchSize;
        }

        bool LearnerSGDEOnOffParallel::checkRefinementNecessary(const std::string &refMonitor, size_t refPeriod,
                                                                size_t totalInstances, double currentValidError,
                                                                double currentTrainError,
                                                                size_t numberOfCompletedRefinements,
                                                                std::shared_ptr<ConvergenceMonitor> &monitor) {
            if (refMonitor == "periodic") {
                // check periodic monitor
                if ((this->offline->getConfig()->decomp_type_ == DBMatDecompChol) &&
                    (totalInstances > 0) && (totalInstances % refPeriod == 0) &&
                    (numberOfCompletedRefinements < this->offline->getConfig()->numRefinements_)) {
                    return true;
                }
            } else if (refMonitor == "convergence") {
                // check convergence monitor
                if (this->validData == nullptr) {
                    //throw base::data_exception(
                    //    "No validation data for checking convergence provided!");
                }
                if ((this->offline->getConfig()->decomp_type_ == DBMatDecompChol) &&
                    (numberOfCompletedRefinements < this->offline->getConfig()->numRefinements_)) {
                    currentValidError = this->getError(*this->validData, *this->validLabels, "Acc");
                    currentTrainError = this->getError(this->trainData, this->trainLabels,
                                                       "Acc");  // if train dataset is large
                    // use a subset for error
                    // evaluation
                    monitor->pushToBuffer(currentValidError, currentTrainError);
                    if (monitor->nextRefCnt > 0) {
                        monitor->nextRefCnt--;
                    }
                    if (monitor->nextRefCnt == 0) {
                        return monitor->checkConvergence();
                    }
                }
            }
            return false;
        }

        size_t
        LearnerSGDEOnOffParallel::handleDataAndZeroBasedRefinement(bool preCompute, MultiGridRefinementFunctor *func,
                                                                   size_t idx, base::Grid *grid,
                                                                   base::GridGenerator &gridGen) const {
            if (preCompute) {
                // precompute the evals (needs to be done once per step, before
                // any refinement is done
                func->preComputeEvaluations();
            }
            func->setGridIndex(idx);
            // perform refinement (zero-crossings-based / data-based)
            size_t gridSizeBeforeRefine = grid->getSize();
            gridGen.refine(*func);
            size_t gridSizeAfterRefine = grid->getSize();
            return gridSizeAfterRefine - gridSizeBeforeRefine;
        }

        size_t LearnerSGDEOnOffParallel::handleSurplusBasedRefinement(DBMatOnlineDE *densEst, base::Grid *grid,
                                                                      base::GridGenerator &gridGen) const {
            //Auxillary Variables
            // required for surplus refinement
            base::DataVector *alphaWork;
            base::DataVector p(trainData.getNcols());

            std::unique_ptr<base::OperationEval> opEval(
                    op_factory::createOperationEval(*grid));
            base::GridStorage &gridStorage = grid->getStorage();
            alphaWork = densEst->getAlpha();
            base::DataVector alphaWeight(alphaWork->getSize());
            // determine surpluses
            for (size_t k = 0; k < gridStorage.getSize(); k++) {
                // sets values of p to the coordinates of the given GridPoint gp
                gridStorage.getPoint(k).getStandardCoordinates(p);
                // multiply k-th alpha with the evaluated function at grind-point
                // k
                alphaWeight[k] = alphaWork->get(k) * opEval->eval(*alphaWork, p);
            }

            // Perform Coarsening (surplus based)
            /*if (coarseCnt < maxCoarseNum) {
              sgpp::base::HashCoarsening coarse_;
              //std::cout << std::endl << "Start coarsening" << std::endl;

              // Coarsening based on surpluses
              sgpp::base::SurplusCoarseningFunctor scf(
                alphaWeight, coarseNumPoints, coarseThreshold);

              //std::cout << "Size before coarsening:" << grid->getSize() <<
            std::endl;
              //int old_size = grid->getSize();
              coarse_.free_coarsen_NFirstOnly(
                grid->getStorage(), scf, alphaWeight, grid->getSize());

              std::cout << "Size after coarsening:" << grid->getSize() <<
            std::endl << std::endl;
              //int new_size = grid->getSize();

              deletedGridPoints.clear();
              deletedGridPoints = coarse_.getDeletedPoints();

              (*refineCoarse)[idx].first = deletedGridPoints;

              coarseCnt++;
            }*/

            // perform refinement (surplus based)
            size_t sizeBeforeRefine = grid->getSize();
            // simple refinement based on surpluses
            base::SurplusRefinementFunctor srf(
                    alphaWeight, offline->getConfig()->ref_noPoints_);
            gridGen.refine(srf);
            size_t sizeAfterRefine = grid->getSize();
            return sizeAfterRefine - sizeBeforeRefine;
        }

        // Train from an entire Batch
        void LearnerSGDEOnOffParallel::train(
                DataBatch *dataBatch,
                bool doCrossValidation,
                std::vector<RefinementResult> *vectorRefinementResults) {
            if (initDone) {
                if (trainData.getNrows() != dataBatch->classLabels->getSize()) {
                    throw sgpp::base::data_exception(
                            "Sizes of train data set and class label vector do not fit!");
                }
                size_t dim = trainData.getNcols();

                // create an empty matrix for every class:
                std::vector<std::pair<sgpp::base::DataMatrix *, double> > trainDataClasses;
                std::map<double, int> classIndices;  // maps class labels to indices

                allocateClassMatrices(dim, trainDataClasses, classIndices);

                // split the data into the different classes:
                splitBatchIntoClasses(dataBatch, dim, trainDataClasses, classIndices);

                // compute density functions
                train(trainDataClasses, doCrossValidation, vectorRefinementResults);

                // delete DataMatrix pointers:
                for (size_t i = 0; i < trainDataClasses.size(); i++) {
                    delete trainDataClasses[i].first;
                }
            }
        }

        void
        LearnerSGDEOnOffParallel::splitBatchIntoClasses(const DataBatch *dataBatch,
                                                        size_t dim,
                                                        const std::vector<std::pair<base::DataMatrix *, double>> &trainDataClasses,
                                                        std::map<double, int> &classIndices) const {
            for (size_t i = 0; i < dataBatch->dataPoints->getNrows(); i++) {
                double classLabel = dataBatch->classLabels->[i];
                base::DataVector vec(dim);
                trainData.getRow(i, vec);
                std::pair<base::DataMatrix *, double> p =
                        trainDataClasses[classIndices[classLabel]];
                p.first->appendRow(vec);
            }
        }

        void LearnerSGDEOnOffParallel::allocateClassMatrices(size_t dim,
                                                             std::vector<std::pair<base::DataMatrix *, double>> &trainDataClasses,
                                                             std::map<double, int> &classIndices) const {
            int index = 0;
            for (size_t i = 0; i < classLabels.getSize(); i++) {
                base::DataMatrix *m = new base::DataMatrix(0, dim);
                std::pair<base::DataMatrix *, double> p(m, classLabels[i]);
                trainDataClasses.push_back(p);
                classIndices.insert(std::pair<double, int>(classLabels[i], index));
                index++;
            }
        }

        //Train from a Batch already split up into its classes
        void LearnerSGDEOnOffParallel::train(
                std::vector<std::pair<sgpp::base::DataMatrix *, double> > &trainDataClasses,
                bool doCrossValidation,
                std::vector<RefinementResult> *vectorRefinementResults) {

            //Calculate the total number of data points
            size_t numberOfDataPoints = 0;
            for (size_t i = 0; i < trainDataClasses.size(); i++) {
                numberOfDataPoints += trainDataClasses[i].first->getSize();
            }

            // Learn from each Class
            for (size_t i = 0; i < trainDataClasses.size(); i++) {
                std::pair<sgpp::base::DataMatrix *, double> p = trainDataClasses[i];

                if ((*p.first).getNrows() > 0) {
                    // update density function for current class
                    (*destFunctions)[i].first->computeDensityFunction(
                            *p.first, true, doCrossValidation, &(*vectorRefinementResults)[i].deletedGridPointsIndexes,
                            (*vectorRefinementResults)[i].numAddedGridPoints);
                    (*vectorRefinementResults)[i].deletedGridPointsIndexes.clear();
                    (*vectorRefinementResults)[i].numAddedGridPoints = 0;

                    if (usePrior) {
                        this->prior[p.second] =
                                ((this->prior[p.second] * static_cast<double>(processedPoints)) +
                                 static_cast<double>(p.first->getSize())) /
                                (static_cast<double>(numberOfDataPoints) +
                                 static_cast<double>(processedPoints));
                    } else {
                        this->prior[p.second] = 1.;
                    }
                }
            }

            this->processedPoints += numberOfDataPoints;
            trained = true;
        }

        double LearnerSGDEOnOffParallel::getAccuracy() {
            sgpp::base::DataVector computedLabels = predict(testData);
            size_t correct = 0;
            size_t correctLabel1 = 0;
            size_t correctLabel2 = 0;
            for (size_t i = 0; i < computedLabels.getSize(); i++) {
                if (computedLabels.get(i) == testLabels.get(i)) {
                    correct++;
                }
                if ((computedLabels.get(i) == -1.0) && (testLabels.get(i) == -1.0)) {
                    correctLabel1++;
                } else if ((computedLabels.get(i) == 1.0) && (testLabels.get(i) == 1.0)) {
                    correctLabel2++;
                }
            }
            // std::cout << "correct label (-1): " << correctLabel1 << std::endl;
            // std::cout << "correct label (1): " << correctLabel2 << std::endl;

            double acc = static_cast<double>(correct) /
                         static_cast<double>(computedLabels.getSize());
            return acc;
        }

        base::DataVector LearnerSGDEOnOffParallel::predict(sgpp::base::DataMatrix &data) {
            base::DataVector result(data.getNrows());

            /*if(not trained) {
              std::cerr << "LearnerSGDEOnOffParallel: Not trained!" << std::endl;
              exit(-1);
            }*/

            for (size_t i = 0; i < data.getNrows(); i++) {
                double max = std::numeric_limits<double>::max() * (-1);
                double max_class = 0;
                // compute the maximum density:
                sgpp::base::DataVector p(data.getNcols());
                data.getRow(i, p);
                for (size_t j = 0; j < destFunctions->size(); j++) {
                    std::pair<DBMatOnlineDE *, double> pair = (*destFunctions)[j];
                    // double density = pair.first->eval(p)*this->prior[pair.second];
                    double density = pair.first->eval(p, true) * this->prior[pair.second];
                    if (density > max) {
                        max = density;
                        max_class = pair.second;
                    }
                }
                if (max_class == 0) {
                    std::cerr << "LearnerSGDEOnOffParallel: Warning: no best class found!"
                              << std::endl;
                }
                result[i] = max_class;
            }
            return result;
        }

        int LearnerSGDEOnOffParallel::predict(sgpp::base::DataVector &p) {
            sgpp::base::DataMatrix ptmp(1, p.getSize());
            ptmp.setRow(0, p);
            sgpp::base::DataVector r = this->predict(ptmp);
            return static_cast<int>(r[0]);
        }

        double LearnerSGDEOnOffParallel::getError(sgpp::base::DataMatrix &data,
                                                  sgpp::base::DataVector &labels,
                                                  std::string errorType) {
            double res = -1.0;

            if (errorType == "Acc") {
                sgpp::base::DataVector computedLabels = predict(data);
                size_t correct = 0;
                for (size_t i = 0; i < computedLabels.getSize(); i++) {
                    if (computedLabels.get(i) == labels.get(i)) {
                        correct++;
                    }
                }

                double acc = static_cast<double>(correct) /
                             static_cast<double>(computedLabels.getSize());
                res = 1.0 - acc;
            }

            return res;
        }

        void LearnerSGDEOnOffParallel::storeResults() {
            sgpp::base::DataVector classesComputed = predict(testData);

            std::ofstream output;
            // write predicted classes to csv file
            output.open("SGDEOnOff_predicted_classes.csv");
            if (output.fail()) {
                std::cout << "failed to create csv file!" << std::endl;
            } else {
                for (size_t i = 0; i < classesComputed.getSize(); i++) {
                    sgpp::base::DataVector x(2);
                    testData.getRow(i, x);
                    output << x[0] << ";" << x[1] << ";" << classesComputed[i] << std::endl;
                }
                output.close();
            }
            // write grids to csv file
            for (size_t i = 0; i < classNumber; i++) {
                DBMatOnlineDE *densEst = (*destFunctions)[i].first;
                sgpp::base::Grid *grid = densEst->getOffline()->getGridPointer();
                output.open("SGDEOnOff_grid_" + std::to_string((*destFunctions)[i].second) +
                            ".csv");
                if (output.fail()) {
                    std::cout << "failed to create csv file!" << std::endl;
                } else {
                    sgpp::base::GridStorage &storage = grid->getStorage();
                    sgpp::base::GridStorage::grid_map_iterator end_iter = storage.end();
                    for (sgpp::base::GridStorage::grid_map_iterator iter = storage.begin();
                         iter != end_iter; iter++) {
                        sgpp::base::DataVector gpCoord(trainData.getNcols());
                        storage.getCoordinates(*(iter->first), gpCoord);
                        for (size_t d = 0; d < gpCoord.getSize(); d++) {
                            if (d < gpCoord.getSize() - 1) {
                                output << gpCoord[d] << ";";
                            } else {
                                output << gpCoord[d] << std::endl;
                            }
                        }
                    }
                    output.close();
                }
            }
            // write density function evaluations to csv file
            double stepSize = 0.01;
            sgpp::base::DataMatrix values(0, 2);
            sgpp::base::DataVector range(101);
            for (size_t i = 0; i < 101; i++) {
                range.set(i, stepSize * (static_cast<double>(i)));
            }
            for (size_t i = 0; i < range.getSize(); i++) {
                for (size_t j = 0; j < range.getSize(); j++) {
                    sgpp::base::DataVector row(2);
                    row.set(1, range.get(i));
                    row.set(0, range.get(j));
                    values.appendRow(row);
                }
            }
            // evaluate each density function at all points from values
            // and write result to csv file
            for (size_t j = 0; j < destFunctions->size(); j++) {
                std::pair<DBMatOnlineDE *, double> pair = (*destFunctions)[j];
                output.open("SGDEOnOff_density_fun_" + std::to_string(pair.second) +
                            "_evals.csv");
                for (size_t i = 0; i < values.getNrows(); i++) {
                    // get next test sample x
                    sgpp::base::DataVector x(2);
                    values.getRow(i, x);
                    double density = pair.first->eval(x, true) * this->prior[pair.second];
                    output << density << ";" << std::endl;
                }
                output.close();
            }
        }

        sgpp::base::DataVector LearnerSGDEOnOffParallel::getDensities(
                sgpp::base::DataVector &point) {
            base::DataVector result(destFunctions->size());
            for (size_t i = 0; i < destFunctions->size(); i++) {
                std::pair<DBMatOnlineDE *, double> pair = (*destFunctions)[i];
                result[i] = pair.first->eval(point);
            }
            return result;
        }

        void LearnerSGDEOnOffParallel::setCrossValidationParameters(
                int lambdaStep, double lambdaStart, double lambdaEnd,
                sgpp::base::DataMatrix *test, sgpp::base::DataMatrix *testRes,
                bool logscale) {
            if (destFunctions != nullptr) {
                for (size_t i = 0; i < destFunctions->size(); i++) {
                    (*destFunctions)[i].first->setCrossValidationParameters(
                            lambdaStep, lambdaStart, lambdaEnd, test, testRes, logscale);
                }
            } else {
                cvSaveLambdaStep = lambdaStep;
                cvSaveLambdaStart = lambdaStart;
                cvSaveLambdaEnd = lambdaEnd;
                cvSaveLogscale = logscale;
                cvSaveTest = test;
                cvSaveTestRes = testRes;
                cvSaved = true;
            }
        }

/*double LearnerSGDEOnOffParallel::getBestLambda() {
  // return online->getBestLambda();
  return 0.; // TODO
}*/

        size_t LearnerSGDEOnOffParallel::getNumClasses() { return this->classNumber; }

        std::vector<std::pair<DBMatOnlineDE *, double> > *
        LearnerSGDEOnOffParallel::getDestFunctions() {
            return destFunctions;
        }

        static const int MPI_MAX_PROCESSOR_NAME_LENGTH = 256;

        void LearnerSGDEOnOffParallel::initMPI() {
            MPI_Init(NULL, NULL);

            // Get World Size
            MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);

            //Get Rank
            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

            //Get Processor Name
            char mpiProcessorName[MPI_MAX_PROCESSOR_NAME_LENGTH];
            int nameLength;
            MPI_Get_processor_name(mpiProcessorName, &nameLength);

            printf("Processor %s (rank %i) has joined MPI pool of size %i\n", mpiProcessorName, world_rank,
                   mpiWorldSize);

            //Setup receiving messages from master/workers
            {
                PendingMPIRequest unicastInputRequest;
                MPI_Packet *mpiPacket = new MPI_Packet;
                unicastInputRequest.buffer = mpiPacket;
                unicastInputRequest.disposeAfterCallback = false;
                unicastInputRequest.callback = [](PendingMPIRequest request) {
                    processIncomingMPICommands(request.buffer);
                    MPI_Irecv(request.buffer, MPI_PACKET_MAX_PAYLOAD_SIZE, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE,
                              MPI_ANY_TAG, MPI_COMM_WORLD, &(request.request));
                };

                MPI_Irecv(mpiPacket, MPI_PACKET_MAX_PAYLOAD_SIZE, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,
                          MPI_COMM_WORLD,
                          &(unicastInputRequest.request));

                pendingMPIRequests.push_back(unicastInputRequest);
                printf("Started listening for unicasts from any sources\n");
            }
            if (!isMaster()) {
                PendingMPIRequest broadcastInputRequest;
                MPI_Packet *mpiPacket = new MPI_Packet;
                broadcastInputRequest.buffer = mpiPacket;
                broadcastInputRequest.disposeAfterCallback = false;
                broadcastInputRequest.callback = [](PendingMPIRequest request) {
                    processIncomingMPICommands(request.buffer);
                    MPI_Ibcast(request.buffer, MPI_PACKET_MAX_PAYLOAD_SIZE, MPI_UNSIGNED_CHAR, MPI_MASTER_RANK,
                               MPI_COMM_WORLD, &(request.request));
                };
                MPI_Ibcast(mpiPacket, MPI_PACKET_MAX_PAYLOAD_SIZE, MPI_UNSIGNED_CHAR, MPI_MASTER_RANK, MPI_COMM_WORLD,
                           &(broadcastInputRequest.request));
                pendingMPIRequests.push_back(broadcastInputRequest);

                printf("Started listening for broadcasts from task master\n");
            }

        }

        void LearnerSGDEOnOffParallel::synchronizeEndOfDataPass() {
            sgpp::parallel::myGlobalMPIComm->Barrier();
        }

        struct RefinementResultNetworkMessage {
            int classIndex;
            int numAddedGridPoints;
            int numDeletedPoints;
            int deletedGridpoints[(MPI_PACKET_MAX_PAYLOAD_SIZE - 3) / sizeof(int)];
        };

        void LearnerSGDEOnOffParallel::sendGridComponentsUpdate(std::vector<RefinementResult> *refinementResults) {
            //TODO Cannot IProbe a Broadcast Message
            sendCommandIDToWorkers(UPDATE_GRID);

            for (size_t classIndex = 0; classIndex < refinementResults->size(); classIndex++) {
                RefinementResult refinementResult = (*refinementResults)[classIndex];

                //Variables to copy block
                size_t numDeletedPoints = 0;
                std::list<size_t>::iterator iterator = refinementResult.deletedGridPointsIndexes.begin();

                for (int deletedPointsIndex = 0; deletedPointsIndex <
                                                 refinementResult.deletedGridPointsIndexes.size(); deletedPointsIndex += numDeletedPoints) {
                    MPI_Packet *mpiPacket = new MPI_Packet;
                    RefinementResultNetworkMessage *networkMessage = (RefinementResultNetworkMessage *) mpiPacket->payload;

                    //First size_t is the class index
                    //Second is the number of added grid points
                    //Third is current number of deletedPoints
                    //More are deleted grid point indices
                    networkMessage->classIndex = classIndex;

                    //The number of added grid points should only be sent the first time
                    networkMessage->numAddedGridPoints =
                            deletedPointsIndex > 0 ? 0 : refinementResult.deletedGridPointsIndexes;

                    numDeletedPoints = std::min(refinementResult.deletedGridPointsIndexes.size() - deletedPointsIndex,
                                                MPI_PACKET_MAX_PAYLOAD_SIZE);

                    //Copy data from list into network buffer
                    for (int bufferIndex = 0; bufferIndex < numDeletedPoints; bufferIndex++) {
                        networkMessage->deletedGridpoints[bufferIndex] = *iterator;
                        iterator++;
                    }

                    //TODO More MPI Packets

                    printf("Sending update for class %zu with %i additions and %zu deletions", classIndex,
                           networkMessage->numAddedGridPoints, numDeletedPoints);
                    PendingMPIRequest pendingMPIRequest;
                    pendingMPIRequest.buffer = mpiPacket;
                    pendingMPIRequests.push_back(pendingMPIRequest);

                    //Send the smallest packet possible
                    MPI_Ibcast(&networkMessage, numDeletedPoints + 3, MPI_UNSIGNED_CHAR, MPI_MASTER_RANK,
                               MPI_COMM_WORLD, &(pendingMPIRequest.request));
                }
            }
        }

        void LearnerSGDEOnOffParallel::processCompletedMPIRequests() {
            MPI_Status mpiStatus;
            int operationCompleted;

            for (std::vector<PendingMPIRequest>::const_iterator pendingMPIRequestIterator = pendingMPIRequests.begin();
                 pendingMPIRequestIterator != pendingMPIRequests.end();
                 pendingMPIRequestIterator++) {
                MPI_Test(&(pendingMPIRequestIterator->request), &operationCompleted, &mpiStatus);
                if (operationCompleted) {
                    //Execute the callback
                    pendingMPIRequestIterator->callback(pendingMPIRequestIterator->buffer);

                    if (pendingMPIRequestIterator->disposeAfterCallback) {
                        //TODO Deleting a void pointer here
                        delete[] pendingMPIRequestIterator->buffer;
                        pendingMPIRequests.erase(pendingMPIRequestIterator);
                    }
                }
            }
        }

        void LearnerSGDEOnOffParallel::waitForMPIRequestsToComplete() {
            for (PendingMPIRequest &pendingMPIRequest : pendingMPIRequests) {
                MPI_Wait(&(pendingMPIRequest.request), MPI_STATUS_IGNORE);
            }

            processCompletedMPIRequests();
        }

        void LearnerSGDEOnOffParallel::receiveGridComponentsUpdate(RefinementResultNetworkMessage *networkMessage) {
            //TODO
            std::list<size_t> deletedGridPointsIndexes;

            for (int deletedGridPoint : networkMessage->deletedGridpoints) {
                deletedGridPointsIndexes.push_back((size_t) deletedGridPoint);
            }

            RefinementResult refinementResult;
            refinementResult.numAddedGridPoints = networkMessage->numAddedGridPoints;
            refinementResult.deletedGridPointsIndexes = deletedGridPointsIndexes;
            updateVariablesAfterRefinement(&refinementResult, networkMessage->classIndex,
                                           (*getDestFunctions())[networkMessage->classIndex].first);
        }


        void LearnerSGDEOnOffParallel::processIncomingMPICommands(MPI_Packet *mpiPacket) {
            switch (mpiPacket->commandID) {
                case UPDATE_GRID:
                    receiveGridComponentsUpdate((RefinementResultNetworkMessage *) (mpiPacket->payload));
                    break;
                case MERGE_GRID:

                    break;
                case ASSIGN_BATCH:
                    break;
                default:
                    printf("Error: MPI unknown command id: %i", mpiPacket->commandID);
                    break;
            }
        }



    }  // namespace datadriven
}  // namespace sgpp

#endif /* USE_GSL */
