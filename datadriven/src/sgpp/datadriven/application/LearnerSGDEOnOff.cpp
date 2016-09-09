// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
//#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>

#include <limits>
#include <ctime>
#include <vector>

LearnerSGDEOnOff::LearnerSGDEOnOff(sgpp::datadriven::DBMatDensityConfiguration& dconf,
                                   sgpp::base::DataMatrix& trainData, sgpp::base::DataVector& trainData_c,
                                   sgpp::base::DataMatrix& testData, sgpp::base::DataVector& testData_c,
                                   /*double* classLabels, */int classNumber, double lambda, bool usePrior, double beta) 
     : destFunctions_(NULL), 
       trained_(false), 
       initDone(false),
       usePrior_(usePrior),  
       beta_(beta) {

  offline = new DBMatOffline(dconf);
  offline->buildMatrix();
  clock_t begin = clock();
  offline->decomposeMatrix();
  clock_t end = clock();
  double elapsed_secs = double(end-begin)/CLOCKS_PER_SEC;
  std::cout << "Decompose matrix: " << elapsed_secs << std::endl;

  // set class labels
  double classLabels[classNumber]; //correct??
  classLabels[0] = -1.0;
  classLabels[1] = 1.0;
  //ToDo: e.g. labels = [-1.0,1.0]
  this->classNumber = classNumber;
  this->classLabels = classLabels;
  
  readOffline(offline); // set offlineObject_ of DBMatOnline class
  //DBMatOnline(offline);

  //If the Cholesky decomposition is chosen declare separate Online-objects for every class
  if (offline->getConfig()->decomp_type_== DBMatDecompChol) {
    destFunctions_ = new vector<std::pair<DBMatOnlineDE*, double> >(classNumber);
    //Every class gets his own online object
    for (int i = 0; i < classNumber; i++) {
      DBMatOnlineDE* densEst = new DBMatOnlineDE(beta_);
      DBMatOffline* offlineRead = new DBMatOffline(*offline);
      densEst->readOffline(offlineRead);
      pair<DBMatOnlineDE*, double> pdest(densEst, classLabels[i]);
      (*destFunctions_)[i] = pdest;
    }
    initDone = true;
  }
  else {
    init();
  }

  this->trainData = &trainData;
  this->trainLabels = &trainData_c;
  this->testData = &testData;
  this->testLabels = &testData_c;
  //this->trainData = std::make_shared<base::DataMatrix>(trainData);
  //this->trainLabels = std::make_shared<base::DataVector>(trainData_c);

  this->cv_saved = false;
  this->processedPoints = 0;

  setLambda(lambda);
}

//void LearnerSGDEOnOff::init(std::map<double, int> entriesPerClass) {
void LearnerSGDEOnOff::init() {
  if (initDone) {
    return;
  }
  //destFunctions_ = new vector<std::pair<DBMatOnlineDE*, double> >(entriesPerClass.size());
  destFunctions_ = new vector<std::pair<DBMatOnlineDE*, double> >(classNumber);
  //int i = 0;
  //for (std::map<double, int>::iterator it = entriesPerClass.begin(); it != entriesPerClass.end(); it++) {
  //std::cout << classNumber << std::endl;
  for (int idx = 0; idx < classNumber; idx++) {
    DBMatOnlineDE* densEst = new DBMatOnlineDE(beta_);
    densEst->readOffline(offline);
    //pair<DBMatOnlineDE*, double> pdest(densEst, (*it).first);
    pair<DBMatOnlineDE*, double> pdest(densEst, classLabels[idx]);
    //(*destFunctions_)[i] = pdest;
    (*destFunctions_)[idx] = pdest;
    //i++;
  }
  if (cv_saved) {
    setCrossValidationParameters(cv_save_lambda_step, cv_save_lambda_start, 
                                 cv_save_lambda_end, cv_save_test, cv_save_test_cc, cv_save_logscale);
  }
  initDone = true;
}

LearnerSGDEOnOff::~LearnerSGDEOnOff() {
  if (destFunctions_ != NULL) {
    for (unsigned int i = 0; i < destFunctions_->size(); i++) {
      /*if ((*destFunctions_)[i].first->getOffline() != NULL) {
        std::cout << "# 1" << std::endl;
        delete (*destFunctions_)[i].first->getOffline();
        std::cout << "# 2" << std::endl;
      }*/
      //std::cout << "# 3" << std::endl;
      delete (*destFunctions_)[i].first;
    }    
    delete destFunctions_;    
  }
}

void LearnerSGDEOnOff::train(size_t batch_size, unsigned int next_cv_step) {
  if(offline->getConfig()->lambda_ < 0) {
    next_cv_step = 1;
  }
  unsigned long long int totalInstances = 0;
  unsigned int next_error_step = 1;
	
  std::pair<sgpp::base::DataMatrix*, sgpp::base::DataVector*> cur_pair;
	
  std::vector<std::pair<std::list<size_t>, unsigned int>>*  RefineCoarse_ = new vector<std::pair<std::list<size_t>, unsigned int>>(classNumber);

  // auxiliary variables
  sgpp::base::DataVector* alpha_work;
  //sgpp::base::GridIndex * gp;
  sgpp::base::DataVector p(trainData->getNcols());

  // dataset passes
  size_t maxIterations = 1;
  size_t numIterations = 0;
  // refinement/coarsening parameters
  size_t ref_period = 25;
  // map number of already performed refinements to grids (classes)
  std::map<int, size_t> refCnts;
  refCnts.insert(std::pair<int, size_t>(-1, 0)); //counter for class -1
  refCnts.insert(std::pair<int, size_t>(1, 0)); //counter for class 1
  // map number of already processed data points to grids (classes)
  std::map<int, size_t> appearances;
  appearances.insert(std::pair<int, size_t>(-1, 0)); //counter for class -1
  appearances.insert(std::pair<int, size_t>(1, 0)); //counter for class 1
  size_t number_coarse = 2;
  double threshold_coarse = 1.0;

  int label;

  std::list<size_t> deletedGridPoints;
  unsigned int newPoints = 0;
	
  std::vector<std::pair<DBMatOnlineDE*, double>>* onlineObjects;

  size_t dim = trainData->getNcols();
  size_t numBatch = trainData->getNrows() / batch_size;

  onlineObjects = getDestFunctions();
  for (int i = 0; i < classNumber; i++) {
    DBMatOnlineDE* densEst = (*onlineObjects)[i].first;
    sgpp::base::Grid* grid = densEst->getOffline()->getGridPointer();
    //size_t  size = grid->getSize();
    std::cout << "initial grid size of grid " << i << " : " << grid->getSize() << std::endl;
  }

  while (numIterations < maxIterations) {
    std::cout << "batches to process: " << numBatch << std::endl;
    size_t cnt = 0;
    for (size_t step = 1; step <= numBatch; step++) {    
      bool do_cv = false;
      if(next_cv_step == step) {
        do_cv = true;
        next_cv_step *= 5;      
      }
      //sgpp::base::DataMatrix batch(0, dim);
      //sgpp::base::DataVector batchLabels(batch_size);
      sgpp::base::DataMatrix* batch = new sgpp::base::DataMatrix(batch_size, dim);
      sgpp::base::DataVector* batchLabels = new sgpp::base::DataVector(batch_size);
      for (size_t j = 0; j < batch_size; j++) {
        sgpp::base::DataVector x(dim);					
        trainData->getRow(j+cnt, x);
        double y = trainLabels->get(j+cnt);
        batch->setRow(j, x);
        batchLabels->set(j, y);
        label = static_cast<int>(y); //this works only for batch size 1 !!!
      }
      appearances.at(label) += 1; //this works only for batch size 1 !!!
      //cur_pair.first = &batch;
      //cur_pair.second = &batchLabels; 
      cur_pair = std::pair<sgpp::base::DataMatrix*, sgpp::base::DataVector*>(batch, batchLabels);
      cnt += batch_size;

      train(*(cur_pair.first), *(cur_pair.second), do_cv, RefineCoarse_);
      totalInstances += (cur_pair.first)->getNrows();
      //std::cout << "Classification before applying adaptivity" << std::endl;
      //doStats(totalInstances, 0, i, next_error_step, false);
    
      //Access DBMatOnlineDE-objects of all classes in order 
      // to apply adaptivity to the specific sparse grids later on
      onlineObjects = getDestFunctions();
      
      //If the Cholesky decomposition is chosen as factorization method refinement 
      //and coarsening methods can be applied
      if ( (offline->getConfig()->decomp_type_== DBMatDecompChol) 
      &&   (appearances.at(label) > 0)
      &&   (appearances.at(label) % ref_period == 0)
      &&   (refCnts.at(label) < offline->getConfig()->numRefinements_) ) {
        std::cout << "refine grid "<< label << std::endl;
        // Bundle grids and surplus vector pointer needed for refinement and
        // evaluation (for zero-crossings refinement, data-based refinement)
        std::vector<sgpp::base::Grid*> grids;
        std::vector<sgpp::base::DataVector*> alphas;
        for (unsigned int i = 0; i < getNumClasses(); i++) {
          DBMatOnlineDE* densEst = (*onlineObjects)[i].first;
          grids.push_back(densEst->getOffline()->getGridPointer());
          alphas.push_back(densEst->getAlpha());
        }
        bool levelPenalize = false;  // Multiplies penalzing term for fine levels
        bool preCompute = true;      // Precomputes and caches evals for zrcr & grid
        sgpp::datadriven::MultiGridRefinementFunctor* func = nullptr;

        // Zero-crossing-based refinement
        sgpp::datadriven::ZeroCrossingRefinementFunctor funcZrcr(grids, alphas,
                                                                 offline->getConfig()->ref_noPoints_,
                                                                 levelPenalize,
                                                                 preCompute);

        // Data-based refinement. Needs a problem dependent coeffA. The values
        // were determined by testing (aim at ~10 % of the training data is
        // to be marked relevant. Cross-validation or similar can/should be employed
        // to determine this value.
        std::vector<double> coeffA;
        coeffA.push_back(1.2); //ripley 1.2
        coeffA.push_back(1.2); //ripley 1.2
        sgpp::datadriven::DataBasedRefinementFunctor funcData(grids, alphas,
                                                              trainData,
                                                              trainLabels,
                                                              offline->getConfig()->ref_noPoints_,
                                                              levelPenalize,
                                                              coeffA);
        func = &funcZrcr;
        //func = &funcData;

        // perform refinement/coarsening for each grid
        //for (unsigned int i = 0; i < getNumClasses(); i++) {
          // perform refinement/coarsening for grid which corresponds to current label
          unsigned int idx;
    	  if (label == -1) {
            idx = 0;
          }
          else if (label == 1) {
            idx = 1;
          }
  	  std::cout << "Refinement and coarsening for class: " << label << std::endl;
          DBMatOnlineDE* densEst = (*onlineObjects)[idx].first;
	  sgpp::base::Grid* grid = densEst->getOffline()->getGridPointer();
          std::unique_ptr<sgpp::base::OperationEval> opEval(op_factory::createOperationEval(*grid));
	  sgpp::base::GridStorage& gridStorage = grid->getStorage();
	  sgpp::base::GridGenerator& gridGen = grid->getGenerator();
		
	  alpha_work = densEst->getAlpha();	
	
	  sgpp::base::DataVector alphaWeight(alpha_work->getSize());

	  //Save grid coordinates
	  /*if(step == 1) {
	    plotGrid(pc, gridStorage, step, i, deletedGridPoints, newPoints, true);
	  }*/
				
	  //Determine surpluses
	  for (unsigned int k = 0; k < gridStorage.getSize(); k++) {
	    //gp = gridStorage.get(k);
            //Sets values of p to the coordinates of the given GridPoint gp
	    //gp->getCoords(p);
            gridStorage.getPoint(k).getStandardCoordinates(p);
	    //Multiply k-th alpha with the evaluated function at grind-point k
	    alphaWeight[k] = alpha_work->get(k)*opEval->eval(*alpha_work, p);
	  } 
          //Performe Coarsening						
	  std::cout << "Size before adaptivity: " << grid->getSize(); 
          /*if (performed_refs > 1) {
	    sgpp::base::HashCoarsening coarse_;
	    //std::cout << std::endl << "Start coarsening" << std::endl;
		
	    //Coarsening based on surpluses
	    sgpp::base::SurplusCoarseningFunctor scf(alphaWeight, number_coarse, threshold_coarse);

	    //std::cout << "Size before coarsening:" << grid->getSize() << std::endl; 
	    //int old_size = grid->getSize();
	    coarse_.free_coarsen_NFirstOnly(grid->getStorage(), scf, alphaWeight, grid->getSize());

	    std::cout << "Size after coarsening:" << grid->getSize() << std::endl << std::endl;
	    //int new_size = grid->getSize();
	  	
	    deletedGridPoints.clear();
	    deletedGridPoints = coarse_.getDeletedPoints();
			
	    (*RefineCoarse_)[idx].first = deletedGridPoints;
          }*/
	  //Perform refinement
	  //std::cout << "Start refinement" << std::endl;
				
	  unsigned int size_before_refine = grid->getSize();

	  //simple refinement based on surpluses
	  //sgpp::base::SurplusRefinementFunctor func(alphaWeight, offline->getConfig()->ref_noPoints_);
	  //gridGen.refine(func);
          
          //zero-crossing / data-based refinement
          if (preCompute) {
            // precompute the evals. Needs to be done once per step, before
            // any refinement is done
            func->preComputeEvaluations();
          }
          gridGen.refine(*func);

          refCnts.at(label) += 1;

	  unsigned int size_after_refine = grid->getSize();
	  std::cout << "  Size after adaptivity: " << grid->getSize() << std::endl;
	  //std::cout << "Size after refinement:" << size_after_refine << std::endl << std::endl;
	  newPoints = size_after_refine - size_before_refine;
	  (*RefineCoarse_)[idx].second = newPoints;
					
	  //plotGrid(gridStorage, step, i, deletedGridPoints, newPoints);

	  //Apply grid changes to the Cholesky factorization	
	  densEst->getOffline()->choleskyModification(newPoints, deletedGridPoints, densEst->getBestLambda());
        //}
      }
      //doStats(totalInstances, 0, step, next_error_step, false);
      delete cur_pair.first;
      delete cur_pair.second;

    }
    numIterations++;
  } //end while
  //trained_ = true;
  //doStats(totalInstances, 0, step, next_error_step, true);
  std::cout << "# Training finished" << std::endl;

  //delete offline;
  delete RefineCoarse_;
}

void LearnerSGDEOnOff::train(sgpp::base::DataMatrix& trainData,
			     sgpp::base::DataVector& classes, bool do_cv, std::vector<std::pair<std::list<size_t>, unsigned int> >*  RefineCoarse_) {
	if (trainData.getNrows() != classes.getSize())
		throw sgpp::base::data_exception(
				"Sizes of train data set and class label vector does not fit together!");
	
	//std::cout << "START SPLITTING UP MATRICES" << std::endl;
	//sgpp::base::SGppStopwatch* myStopwatch = new sgpp::base::SGppStopwatch();
	//	myStopwatch->start();
	
	int dim = trainData.getNcols();

	//Compute all occurring class labels and how many data points exist per label:
	std::map<double, int> entriesPerClass;
        //if(!initDone) {
	//  for (size_t i = 0; i < classes.getSize(); i++) {
	//    double classNum = classes.get(i);
	    /*if (entriesPerClass.find(classNum) != entriesPerClass.end()) {
	      entriesPerClass.insert(std::pair<double, int>(classNum, 1));
	      std::cout << "new class: " << classNum << std::endl;
	    } 
            else {
	      entriesPerClass[classNum]++; //TODO get this running
	    }*/
	//    entriesPerClass[classNum]++;
	//  }
	//  init(entriesPerClass);
	//}
        //else {
	  for (size_t i = 0; i < destFunctions_->size(); i++) {
	    entriesPerClass[(*destFunctions_)[i].second] = 1; // We only need the key set, the value is not important here
	  }
        //}	
	//Create an empty matrix for every class:
	std::vector<std::pair<sgpp::base::DataMatrix*, double> > trainDataClasses;
	std::map<double, int> class_indeces; //Maps class numbers to indices

	std::map<double, int>::iterator it;
	int index = 0;
	for (it = entriesPerClass.begin(); it != entriesPerClass.end(); it++) {
		sgpp::base::DataMatrix* m = new sgpp::base::DataMatrix(0/*(*it).second*/,
				dim);
		pair<sgpp::base::DataMatrix*, double> p(m, (*it).first);
		trainDataClasses.push_back(p);
		class_indeces[(*it).first] = index;
		index++;
	}
	//Split the data into the different classes:
	for (size_t i = 0; i < trainData.getNrows(); i++) {
		double classLabel = classes[i];
		sgpp::base::DataVector vec(dim);
		trainData.getRow(i, vec);
		pair<sgpp::base::DataMatrix*, double> p =
				trainDataClasses[class_indeces[classLabel]];
		p.first->appendRow(vec);
	}
	
	//std::cout << "Time: " << myStopwatch->stop() << std::endl;
	train(trainDataClasses, do_cv, RefineCoarse_);

	//delete DataMatrix pointers:
	for (size_t i = 0; i < trainDataClasses.size(); i++) {
		delete trainDataClasses[i].first;
	}
}

void LearnerSGDEOnOff::train( std::vector<std::pair<sgpp::base::DataMatrix*, double> >& trainDataClasses, bool do_cv, std::vector<std::pair<std::list<size_t>, unsigned int> >*  RefineCoarse_) {

	time_t startTime = time(NULL);

	unsigned int numberOfDataPoints = 0;
	for(unsigned int i = 0; i < trainDataClasses.size(); i++)
		numberOfDataPoints += trainDataClasses[i].first->getSize();

	for (unsigned int i = 0; i < trainDataClasses.size(); i++) {
		pair<sgpp::base::DataMatrix*, double> p = trainDataClasses[i];

		//Create OnlineDE object:
		(*destFunctions_)[i].first->computeDensityFunction(*p.first, true, do_cv, &(*RefineCoarse_)[i].first, (*RefineCoarse_)[i].second);

		//(*destFunctions_)[i].first->normalize(1000);
		if(usePrior_){
		  this->prior[p.second] = ((this->prior[p.second] * processedPoints) + p.first->getSize())/((double)numberOfDataPoints + processedPoints);
		}else
		  this->prior[p.second] = 1.;
	}

	this->processedPoints += numberOfDataPoints;

	this->traintime = difftime(time(NULL), startTime);
	//std::cout << "Time to train: " << this->traintime << std::endl;
	//std::cout << "Time: " << myStopwatch->stop() << std::endl;

	trained_ = true;
}

double LearnerSGDEOnOff::getAccuracy() {
  sgpp::base::DataVector computedLabels = predict(testData);
  int correct = 0;
  for (int i = 0; i < computedLabels.getSize(); i++) {
    if (computedLabels.get(i) == testLabels->get(i)) {
      correct++;
    }
  }
  double acc = static_cast<double>(correct) / static_cast<double>(computedLabels.getSize());
  return acc;
}

sgpp::base::DataVector LearnerSGDEOnOff::predict(sgpp::base::DataMatrix* data) {
	sgpp::base::DataVector result(data->getNrows());

	if(not trained_) {
	  std::cerr << "LearnerSGDEOnOff: Not trained!" << std::endl;
	  exit(-1);
	}

	for (unsigned int i = 0; i < data->getNrows(); i++) {
		double max = numeric_limits<double>::max()*(-1);
		double max_class = 0;
		//Compute the maximum density:
		sgpp::base::DataVector p(data->getNcols());
		data->getRow(i, p);
		for (unsigned int j = 0; j < destFunctions_->size(); j++) {
			pair<DBMatOnlineDE*, double> pair = (*destFunctions_)[j];
			double density = pair.first->eval(p)*this->prior[pair.second];
			if (density > max) {
				max = density;
				max_class = pair.second;
			}
		}
		//if(max_class == -1)
		//	std::cerr << "LearnerSGDEOnOff: Warning: no best class found!" << std::endl;
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

sgpp::base::DataVector LearnerSGDEOnOff::getDensities(sgpp::base::DataVector& point) {

	sgpp::base::DataVector result (destFunctions_->size());
	for (unsigned int i = 0; i < destFunctions_->size(); i++) {
		pair<DBMatOnlineDE*, double> pair = (*destFunctions_)[i];
		result[i] = pair.first->eval(point);
	}
	return result;
}

void LearnerSGDEOnOff::setCrossValidationParameters(int lambda_step, double lambda_start, double lambda_end,
                                               sgpp::base::DataMatrix *test, sgpp::base::DataMatrix *test_cc, 
                                               bool logscale) {
	if (destFunctions_ != NULL) {
		for (size_t i = 0; i < destFunctions_->size(); i++) {
			(*destFunctions_)[i].first->setCrossValidationParameters(lambda_step, lambda_start, lambda_end, test, test_cc, logscale);
		}
	}
	else {
		cv_save_lambda_step = lambda_step;
		cv_save_lambda_start = lambda_start;
		cv_save_lambda_end = lambda_end;
		cv_save_logscale = logscale;
		cv_save_test = test;
		cv_save_test_cc = test_cc;
		cv_saved = true;
	}
}

double LearnerSGDEOnOff::getBestLambda() {
  //return online->getBestLambda();
  return 0.; // TODO
}

unsigned int LearnerSGDEOnOff::getNumClasses() {
  return this->classNumber;
}

/*void LearnerSGDEOnOff::normalize(size_t samples) {
	for (size_t i = 0; i < destFunctions_->size(); i++) {
		std::cout << std::endl << (*destFunctions_)[i].first->normalize(samples) << std::endl;
	}
}*/

std::vector<std::pair<DBMatOnlineDE*, double> >* LearnerSGDEOnOff::getDestFunctions(){
  return destFunctions_;
}
