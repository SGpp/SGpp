// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/algorithm/ConvergenceMonitor.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>

#include <limits>
#include <ctime>
#include <vector>
#include <cmath>

namespace sgpp {
namespace datadriven {

LearnerSGDEOnOff::LearnerSGDEOnOff(sgpp::datadriven::DBMatDensityConfiguration& dconf,
                                   sgpp::base::DataMatrix& trainData, 
                                   sgpp::base::DataVector& trainDataLabels,
                                   sgpp::base::DataMatrix& testData, 
                                   sgpp::base::DataVector& testDataLabels,
                                   sgpp::base::DataMatrix* validData, 
                                   sgpp::base::DataVector* validDataLabels,
                                   sgpp::base::DataVector& classLabels, unsigned int classNumber, 
                                   bool usePrior, double beta, double lambda) 
     : trainData(&trainData),
       trainLabels(&trainDataLabels),
       testData(&testData),
       testLabels(&testDataLabels),
       validData(validData),
       validLabels(validDataLabels),
       classLabels(classLabels),
       classNumber(classNumber), 
       trained(false), 
       initDone(false),
       usePrior(usePrior),  
       beta(beta),
       destFunctions(nullptr),
       error(-1.0) {

  offline = new DBMatOffline(dconf);
  offline->buildMatrix();
  //clock_t begin = clock();
  offline->decomposeMatrix();
  //clock_t end = clock();
  //double elapsed_secs = double(end-begin)/CLOCKS_PER_SEC;
  //std::cout << "#Decompose matrix: " << elapsed_secs << std::endl;*/
  
  readOffline(offline); // set offlineObject_ of DBMatOnline class

  // if the Cholesky decomposition is chosen declare separate Online-objects for every class
  if (offline->getConfig()->decomp_type_== DBMatDecompChol) {
    destFunctions = new vector<std::pair<DBMatOnlineDE*, double>>(classNumber);
    // every class gets his own online object
    for (unsigned int i = 0; i < classNumber; i++) {
      DBMatOnlineDE* densEst = new DBMatOnlineDE(beta);
      DBMatOffline* offlineRead = new DBMatOffline(*offline);
      densEst->readOffline(offlineRead);
      pair<DBMatOnlineDE*, double> pdest(densEst, classLabels[i]);
      (*destFunctions)[i] = pdest;
    }
    initDone = true;
  }
  else {
    init();
  }

  cvSaved = false;
  processedPoints = 0;

  for (unsigned int i = 0; i < classNumber; i++) {
    prior.insert(std::pair<double, double>(classLabels[i], 0.0));
  } 

  setLambda(lambda);
}

void LearnerSGDEOnOff::init() {
  if (initDone) {
    return;
  }
  destFunctions = new vector<std::pair<DBMatOnlineDE*, double>>(classNumber);
  for (unsigned int idx = 0; idx < classNumber; idx++) {
    DBMatOnlineDE* densEst = new DBMatOnlineDE(beta);
    densEst->readOffline(offline);
    pair<DBMatOnlineDE*, double> pdest(densEst, classLabels[idx]);
    (*destFunctions)[idx] = pdest;
  }
  if (cvSaved) {
    setCrossValidationParameters(cvSaveLambdaStep, cvSaveLambdaStart, 
                                 cvSaveLambdaEnd, cvSaveTest, 
                                 cvSaveTestRes, cvSaveLogscale);
  }
  initDone = true;
}

LearnerSGDEOnOff::~LearnerSGDEOnOff() {
  if (destFunctions != nullptr) {
    for (unsigned int i = 0; i < destFunctions->size(); i++) {
      delete (*destFunctions)[i].first;
    }    
    delete destFunctions;    
  }
}

void LearnerSGDEOnOff::train(size_t batchSize, size_t maxDataPasses, 
                             string refType, string refMonitor, 
                             size_t refPeriod, double accDeclineThreshold, 
                             size_t accDeclineBufferSize, size_t minRefInterval,
                             bool enableCv, unsigned int nextCvStep) {
  // counts total number of processed data points
  size_t totalInstances = 0;
  // pointer to the next batch (data points + class labels) to be processed
  std::pair<sgpp::base::DataMatrix*, sgpp::base::DataVector*> curPair;
  // contains list of removed grid points and number of added grid points
  // (is updated in each refinement/coarsening step)	
  std::vector<std::pair<std::list<size_t>, unsigned int>>*  refineCoarse = 
    new vector<std::pair<std::list<size_t>, unsigned int>>(classNumber);

  // auxiliary variables
  sgpp::base::DataVector* alphaWork; // required for surplus refinement
  sgpp::base::DataVector p(trainData->getNcols());

  // initialize counter for dataset passes 
  size_t cntDataPasses = 0;

  // initialize refinement variables 
  double currentValidError = 0.0;
  double currentTrainError = 0.0;
  // create convergence monitor object
  std::shared_ptr<ConvergenceMonitor> monitor(new ConvergenceMonitor(
    accDeclineThreshold, accDeclineBufferSize, minRefInterval));
  bool doRefine = false; // is set to 'true' by refinement monitor
  // counts number of performed refinements
  size_t refCnt = 0;

  // coarsening
  //size_t coarseCnt = 0;
  //size_t maxCoarseNum = 5;
  //size_t coarsePeriod = 50;
  //size_t coarseNumPoints = 1;
  //double coarseThreshold = 1.0;

  std::list<size_t> deletedGridPoints;
  unsigned int newPoints = 0;
	
  std::vector<std::pair<DBMatOnlineDE*, double>>* onlineObjects;

  size_t dim = trainData->getNcols();
  // determine number of batches to process
  size_t numBatch = trainData->getNrows() / batchSize;

  // print initial grid size
  onlineObjects = getDestFunctions();
  for (unsigned int i = 0; i < classNumber; i++) {
    DBMatOnlineDE* densEst = (*onlineObjects)[i].first;
    sgpp::base::Grid* grid = densEst->getOffline()->getGridPointer();
    std::cout << "#Initial grid size of grid " << i << " : " << grid->getSize() << std::endl;
  }

  // auxiliary variable for accuracy (error) measurement 
  double acc = 0.0;
  //acc = getAccuracy();
  //avgErrors.append(1.0 - acc);  


  // main loop which performs the training process
  while (cntDataPasses < maxDataPasses) {
    std::cout << "#batch-size: " << batchSize << std::endl;
    std::cout << "#batches to process: " << numBatch << std::endl;
    size_t cnt = 0; // data point counter - determines offset when selecting next batch
    // iterate over total number of batches 
    for (size_t step = 1; step <= numBatch; step++) {  
      // check if cross-validation should be performed
      bool doCv = false;
      if (enableCv) {  
        if(nextCvStep == step) {
          doCv = true;
          nextCvStep *= 5;  
        }         
      }
      // assemble next batch
      sgpp::base::DataMatrix* batch = new sgpp::base::DataMatrix(batchSize, dim);
      sgpp::base::DataVector* batchLabels = new sgpp::base::DataVector(batchSize);
      for (size_t j = 0; j < batchSize; j++) {
        sgpp::base::DataVector x(dim);					
        trainData->getRow(j+cnt, x);
        double y = trainLabels->get(j+cnt);
        batch->setRow(j, x);
        batchLabels->set(j, y);
      }
      curPair = std::pair<sgpp::base::DataMatrix*, sgpp::base::DataVector*>(batch, batchLabels);
      cnt += batchSize;
      
      // train the model with current batch
      train(*(curPair.first), *(curPair.second), doCv, refineCoarse);

      totalInstances += (curPair.first)->getNrows();
    
      // access DBMatOnlineDE-objects of all classes in order 
      // to apply adaptivity to the specific sparse grids later on
      onlineObjects = getDestFunctions();

      // check if refinement should be performed
      if (refMonitor == "periodic") {
        // check periodic monitor
        if ( (offline->getConfig()->decomp_type_== DBMatDecompChol)
        &&   (totalInstances > 0)
        &&   (totalInstances % refPeriod == 0)
        &&   (refCnt < offline->getConfig()->numRefinements_) ) {
          doRefine = true;
        }
      }
      else if (refMonitor == "convergence") {
        // check convergence monitor
        if (validData == nullptr) {
          throw base::data_exception(
	    "No validation data for checking convergence provided!");
        }
        if ( (offline->getConfig()->decomp_type_== DBMatDecompChol)  
        &&   (refCnt < offline->getConfig()->numRefinements_) ) {
          currentValidError = getError(validData, validLabels, "Acc");
          currentTrainError = getError(trainData, trainLabels, "Acc"); // if train dataset is large
                                                                       // use a subset for error 
                                                                       // evaluation
          monitor->pushToBuffer(currentValidError,currentTrainError);
          if (monitor->nextRefCnt > 0) {
            monitor->nextRefCnt--;
          }          
          if (monitor->nextRefCnt == 0) {
            doRefine = monitor->checkConvergence();
          } 
        }
      }
      // if the Cholesky decomposition is chosen as factorization method refinement 
      // and coarsening methods can be applied
      if (doRefine) {
        //acc = getAccuracy();
        //avgErrors.append(1.0 - acc);
        std::cout << "refinement at iteration: "<< totalInstances << std::endl;
        // bundle grids and surplus vector pointer needed for refinement
        // (for zero-crossings refinement, data-based refinement)
        std::vector<sgpp::base::Grid*> grids;
        std::vector<sgpp::base::DataVector*> alphas;
        for (unsigned int i = 0; i < getNumClasses(); i++) {
          DBMatOnlineDE* densEst = (*onlineObjects)[i].first;
          grids.push_back(densEst->getOffline()->getGridPointer());
          alphas.push_back(densEst->getAlpha());
        }
        bool levelPenalize = false;  // Multiplies penalzing term for fine levels
        bool preCompute = true;      // Precomputes and caches evals for zrcr 
        sgpp::datadriven::MultiGridRefinementFunctor* func = nullptr;

        
        // Zero-crossing-based refinement
        sgpp::datadriven::ZeroCrossingRefinementFunctor funcZrcr =  
          *(new sgpp::datadriven::ZeroCrossingRefinementFunctor(grids, alphas,
                                                                offline->getConfig()->ref_noPoints_,
                                                                levelPenalize,
                                                                preCompute));
            
        // Data-based refinement. Needs a problem dependent coeffA. The values
        // can be determined by testing (aim at ~10 % of the training data is
        // to be marked relevant). Cross-validation or similar can/should be employed
        // to determine this value.
        std::vector<double> coeffA;
        coeffA.push_back(1.2); //ripley 1.2
        coeffA.push_back(1.2); //ripley 1.2
        sgpp::datadriven::DataBasedRefinementFunctor funcData = 
          *(new sgpp::datadriven::DataBasedRefinementFunctor(grids, alphas,
                                                             trainData,
                                                             trainLabels,
                                                             offline->getConfig()->ref_noPoints_,
                                                             levelPenalize,
                                                             coeffA));  
        if (refType == "zero") {
          func = &funcZrcr;
        }
        else if (refType == "data") {
          func = &funcData;
        }
        
        // perform refinement/coarsening for each grid
        unsigned int sizeBeforeRefine;
        unsigned int sizeAfterRefine;
        for (size_t idx = 0; idx < getNumClasses(); idx++) {
          //perform refinement/coarsening for grid which corresponds to current index
  	  std::cout << "Refinement and coarsening for class: " << idx << std::endl;
          DBMatOnlineDE* densEst = (*onlineObjects)[idx].first;
	  sgpp::base::Grid* grid = densEst->getOffline()->getGridPointer();
          std::cout << "Size before adaptivity: " << grid->getSize() << std::endl;
          
          sgpp::base::GridGenerator& gridGen = grid->getGenerator();

          if (refType == "surplus") {
            std::unique_ptr<sgpp::base::OperationEval> opEval(
              op_factory::createOperationEval(*grid));
	    sgpp::base::GridStorage& gridStorage = grid->getStorage();	  	
	    alphaWork = densEst->getAlpha();		
	    sgpp::base::DataVector alphaWeight(alphaWork->getSize());			
	    // determine surpluses
	    for (unsigned int k = 0; k < gridStorage.getSize(); k++) {
              // sets values of p to the coordinates of the given GridPoint gp
              gridStorage.getPoint(k).getStandardCoordinates(p);
	      // multiply k-th alpha with the evaluated function at grind-point k
	      alphaWeight[k] = alphaWork->get(k)*opEval->eval(*alphaWork, p);
	    }

            // Perform Coarsening (surplus based)
            /*if (coarseCnt < maxCoarseNum) {
	      sgpp::base::HashCoarsening coarse_;
	      //std::cout << std::endl << "Start coarsening" << std::endl;
		
	      // Coarsening based on surpluses
	      sgpp::base::SurplusCoarseningFunctor scf(
                alphaWeight, coarseNumPoints, coarseThreshold);

	      //std::cout << "Size before coarsening:" << grid->getSize() << std::endl; 
	      //int old_size = grid->getSize();
	      coarse_.free_coarsen_NFirstOnly(
                grid->getStorage(), scf, alphaWeight, grid->getSize());

	      std::cout << "Size after coarsening:" << grid->getSize() << std::endl << std::endl;
	      //int new_size = grid->getSize();
	  	
	      deletedGridPoints.clear();
	      deletedGridPoints = coarse_.getDeletedPoints();
			
	      (*refineCoarse)[idx].first = deletedGridPoints;

              coarseCnt++;
            }*/ 

	    // perform refinement (surplus based)				
	    sizeBeforeRefine = grid->getSize();
	    // simple refinement based on surpluses
	    sgpp::base::SurplusRefinementFunctor srf(alphaWeight, 
                                                     offline->getConfig()->ref_noPoints_);
	    gridGen.refine(srf);
            sizeAfterRefine = grid->getSize();
          }
          else if ( (refType == "data") || (refType == "zero") ) {
            if (preCompute) {
              // precompute the evals (needs to be done once per step, before
              // any refinement is done
              func->preComputeEvaluations();
            }
            func->setGridIndex(idx);
            // perform refinement (zero-crossings-based / data-based)
            sizeBeforeRefine = grid->getSize();
            gridGen.refine(*func);
            sizeAfterRefine = grid->getSize();
          }
	  
	  std::cout << "grid size after adaptivity: " << grid->getSize() << std::endl;

	  newPoints = sizeAfterRefine - sizeBeforeRefine;
	  (*refineCoarse)[idx].second = newPoints;

	  // apply grid changes to the Cholesky factorization	
	  densEst->getOffline()->choleskyModification(
            newPoints, deletedGridPoints, densEst->getBestLambda());
          // update alpha vector
          densEst->updateAlpha(&(*refineCoarse)[idx].first, (*refineCoarse)[idx].second);
        }
        refCnt += 1;
        doRefine = false;
        if (refMonitor == "convergence") {
          monitor->nextRefCnt = monitor->minRefInterval;
        }
      }
      delete curPair.first;
      delete curPair.second;
     
      // save current error
      /*if (totalInstances % 10 == 0) {
        acc = getAccuracy();
        avgErrors.append(1.0 - acc);  
      }*/

    }
    cntDataPasses++;
    processedPoints = 0;

  } // end while

  std::cout << "#Training finished" << std::endl;

  error = 1.0 - getAccuracy();

  //delete offline;
  delete refineCoarse;
}

void LearnerSGDEOnOff::train(sgpp::base::DataMatrix& trainData,
			     sgpp::base::DataVector& trainClasses, bool doCv, 
                             std::vector<std::pair<std::list<size_t>,unsigned int>>* refineCoarse) {
  if (initDone) {
    if (trainData.getNrows() != trainClasses.getSize()) {
      throw sgpp::base::data_exception(
        "Sizes of train data set and class label vector do not fit!");
    }
    size_t dim = trainData.getNcols();

    // create an empty matrix for every class:
    std::vector<std::pair<sgpp::base::DataMatrix*, double>> trainDataClasses;
    std::map<double, int> classIndices; // maps class labels to indices
    int index = 0;
    for (size_t i = 0; i < classLabels.getSize(); i++) {
      sgpp::base::DataMatrix* m = new sgpp::base::DataMatrix(0, dim);
      pair<sgpp::base::DataMatrix*, double> p(m, classLabels[i]);
      trainDataClasses.push_back(p);
      classIndices.insert(std::pair<double, int>(classLabels[i], index));
      index++;
    }
    // split the data into the different classes:
    for (size_t i = 0; i < trainData.getNrows(); i++) {
      double classLabel = trainClasses[i];
      sgpp::base::DataVector vec(dim);
      trainData.getRow(i, vec);
      pair<sgpp::base::DataMatrix*, double> p =
        trainDataClasses[classIndices[classLabel]];
      p.first->appendRow(vec);
    }
    // compute density functions
    train(trainDataClasses, doCv, refineCoarse);

    // delete DataMatrix pointers:
    for (size_t i = 0; i < trainDataClasses.size(); i++) {
      delete trainDataClasses[i].first;
    }
  }
}

void LearnerSGDEOnOff::train(
  std::vector<std::pair<sgpp::base::DataMatrix*, double>>& trainDataClasses, 
  bool doCv, 
  std::vector<std::pair<std::list<size_t>, unsigned int>>*  refineCoarse) {

  unsigned int numberOfDataPoints = 0;
  for (unsigned int i = 0; i < trainDataClasses.size(); i++) {
    numberOfDataPoints += trainDataClasses[i].first->getSize();
  }
  for (unsigned int i = 0; i < trainDataClasses.size(); i++) {
    pair<sgpp::base::DataMatrix*, double> p = trainDataClasses[i];

    if ((*p.first).getNrows() > 0) {
      // update density function for current class
      (*destFunctions)[i].first->computeDensityFunction(
        *p.first, true, doCv, &(*refineCoarse)[i].first, (*refineCoarse)[i].second);
      (*refineCoarse)[i].first.clear(); 
      (*refineCoarse)[i].second = 0;

      if (usePrior) {
        this->prior[p.second] = ((this->prior[p.second] * processedPoints) + p.first->getSize()) /
          ((double)numberOfDataPoints + processedPoints);
      } 
      else {
        this->prior[p.second] = 1.;
      }
    }
  }

  this->processedPoints += numberOfDataPoints;
  trained = true;
}

double LearnerSGDEOnOff::getAccuracy() {
  sgpp::base::DataVector computedLabels = predict(testData);
  int correct = 0;
  int correctLabel1 = 0;
  int correctLabel2 = 0;
  for (int i = 0; i < computedLabels.getSize(); i++) {
    if (computedLabels.get(i) == testLabels->get(i)) {
      correct++;
    }
    if ((computedLabels.get(i) == -1.0) && (testLabels->get(i) == -1.0)) {
      correctLabel1++;
    }
    else if ((computedLabels.get(i) == 1.0) && (testLabels->get(i) == 1.0)) {
      correctLabel2++;
    }
  }
  //std::cout << "correct label (-1): " << correctLabel1 << std::endl;
  //std::cout << "correct label (1): " << correctLabel2 << std::endl;

  double acc = static_cast<double>(correct) / static_cast<double>(computedLabels.getSize());
  return acc;
}

base::DataVector LearnerSGDEOnOff::predict(sgpp::base::DataMatrix* data) {
  base::DataVector result(data->getNrows());

  /*if(not trained) {
    std::cerr << "LearnerSGDEOnOff: Not trained!" << std::endl;
    exit(-1);
  }*/

  for (unsigned int i = 0; i < data->getNrows(); i++) {
    double max = numeric_limits<double>::max()*(-1);
    double max_class = 0;
    // compute the maximum density:
    sgpp::base::DataVector p(data->getNcols());
    data->getRow(i, p);
    for (unsigned int j = 0; j < destFunctions->size(); j++) {
      pair<DBMatOnlineDE*, double> pair = (*destFunctions)[j];
      //double density = pair.first->eval(p)*this->prior[pair.second];
      double density = pair.first->eval(p, true)*this->prior[pair.second];
      if (density > max) {
        max = density;
        max_class = pair.second;
      }
    }
    if(max_class == 0) {
      std::cerr << "LearnerSGDEOnOff: Warning: no best class found!" << std::endl;
    }
    result[i] = max_class;
  }
  return result;
}

int LearnerSGDEOnOff::predict(sgpp::base::DataVector &p) {
  sgpp::base::DataMatrix ptmp(1, p.getSize());
  ptmp.setRow(0, p);
  sgpp::base::DataVector r = this->predict(&ptmp);
  return static_cast<int>(r[0]);
}

double LearnerSGDEOnOff::getError(sgpp::base::DataMatrix* data, sgpp::base::DataVector* labels, 
                                  std::string errorType) {
  double res = -1.0;

  if (errorType == "Acc") {
    sgpp::base::DataVector computedLabels = predict(data);
    int correct = 0;
    for (int i = 0; i < computedLabels.getSize(); i++) {
      if (computedLabels.get(i) == labels->get(i)) {
        correct++;
      }
    }

    double acc = static_cast<double>(correct) / static_cast<double>(computedLabels.getSize());
    res = 1.0 - acc;
  }
  
  return res;
}

void LearnerSGDEOnOff::storeResults() {

  sgpp::base::DataVector classesComputed = predict(testData);

  std::ofstream output;
  // write predicted classes to csv file
  output.open("SGDEOnOff_predicted_classes.csv");
  if (output.fail()) {
    std::cout << "failed to create csv file!" << std::endl;  
  }
  else {
    for (size_t i = 0; i < classesComputed.getSize(); i++) {
      sgpp::base::DataVector x(2);					
      testData->getRow(i, x);
      output << x[0] << ";" << x[1] << ";" << classesComputed[i] << std::endl;
    }
    output.close();
  }
  // write grids to csv file
  for (unsigned int i = 0; i < classNumber; i++) {
    DBMatOnlineDE* densEst = (*destFunctions)[i].first;
    sgpp::base::Grid* grid = densEst->getOffline()->getGridPointer();  
    output.open("SGDEOnOff_grid_"+std::to_string((*destFunctions)[i].second)+".csv");
    if (output.fail()) {
      std::cout << "failed to create csv file!" << std::endl;  
    }
    else {
      sgpp::base::GridStorage& storage = grid->getStorage();
      sgpp::base::GridStorage::grid_map_iterator end_iter = storage.end();
      for (sgpp::base::GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter; iter++) { 
        sgpp::base::DataVector gpCoord(trainData->getNcols());
        storage.getCoordinates(*(iter->first), gpCoord);
        for (size_t d = 0; d < gpCoord.getSize(); d++) {
          if (d < gpCoord.getSize()-1) {
            output << gpCoord[d] << ";";
          }
          else {
            output << gpCoord[d] << std::endl;
          }
        }
      }
      output.close();
    }
  }
  // write density function evaluations to csv file
  double stepSize = 0.01;
  sgpp::base::DataMatrix values(0,2);
  sgpp::base::DataVector range(101);
  for (size_t i = 0; i < 101; i++) {
    range.set(i, stepSize*(static_cast<double>(i)));
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
  for (unsigned int j = 0; j < destFunctions->size(); j++) {
    pair<DBMatOnlineDE*, double> pair = (*destFunctions)[j];
    output.open("SGDEOnOff_density_fun_"+std::to_string(pair.second)+"_evals.csv");
    for (size_t i = 0; i < values.getNrows(); i++) {
      // get next test sample x 
      sgpp::base::DataVector x(2);					
      values.getRow(i, x);
      double density = pair.first->eval(x, true)*this->prior[pair.second];
      output << density << ";" << std::endl;
    }
    output.close();
  }

}

sgpp::base::DataVector LearnerSGDEOnOff::getDensities(sgpp::base::DataVector& point) {
  base::DataVector result(destFunctions->size());
  for (size_t i = 0; i < destFunctions->size(); i++) {
    pair<DBMatOnlineDE*, double> pair = (*destFunctions)[i];
    result[i] = pair.first->eval(point);
  }
  return result;
}

void LearnerSGDEOnOff::setCrossValidationParameters(int lambdaStep, double lambdaStart, 
                                                    double lambdaEnd, sgpp::base::DataMatrix* test, 
                                                    sgpp::base::DataMatrix* testRes, 
                                                    bool logscale) {
  if (destFunctions != nullptr) {
    for (size_t i = 0; i < destFunctions->size(); i++) {
      (*destFunctions)[i].first->setCrossValidationParameters(
        lambdaStep, lambdaStart, lambdaEnd, test, testRes, logscale);
    }
  }
  else {
    cvSaveLambdaStep = lambdaStep;
    cvSaveLambdaStart = lambdaStart;
    cvSaveLambdaEnd = lambdaEnd;
    cvSaveLogscale = logscale;
    cvSaveTest = test;
    cvSaveTestRes = testRes;
    cvSaved = true;
  }
}

/*double LearnerSGDEOnOff::getBestLambda() {
  //return online->getBestLambda();
  return 0.; // TODO
}*/

unsigned int LearnerSGDEOnOff::getNumClasses() {
  return this->classNumber;
}

std::vector<std::pair<DBMatOnlineDE*, double>>* LearnerSGDEOnOff::getDestFunctions(){
  return destFunctions;
}

}  // namespace datadriven
}  // namespace sgpp
