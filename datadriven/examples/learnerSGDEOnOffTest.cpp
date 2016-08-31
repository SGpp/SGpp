// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>

#include <ctime>
#include <string>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

int main(int argc, char** argv) {
  std::string filename = "../tests/data/ripley/ripley_train_1.arff";
  //std::string filename = "../tests/data/ripleyGarcke.train.arff";
  //std::string filename = "../tests/data/banana/banana_train_1.arff";
  //std::string filename = "../tests/data/banana.arff";
  // load training samples
  std::cout << "# loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset trainDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
  sgpp::base::DataMatrix& trainData = trainDataset.getData();

  // normalize data
  trainData.normalizeDimension(0);
  trainData.normalizeDimension(1);

  // extract train classes
  sgpp::base::DataVector& trainLabels = trainDataset.getTargets();

  filename = "../tests/data/ripley/ripley_test.arff";
  //filename = "../tests/data/ripleyGarcke.test.arff";
  //filename = "../tests/data/banana/banana_train_0.arff";
  //filename = "../tests/data/banana.arff";
  // load test samples
  std::cout << "# loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset testDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
  sgpp::base::DataMatrix& testData = testDataset.getData();

  // normalize data
  testData.normalizeDimension(0);
  testData.normalizeDimension(1);

  // extract test classes
  sgpp::base::DataVector& testLabels = testDataset.getTargets();

  // set number of classes
  int classNum = 2;
 
  // set class labels
  /*double labels[classNum]; //correct??
  labels[0] = -1.0;
  labels[1] = 1.0;*/
  //ToDo: e.g. labels = [-1.0,1.0]
  

  // configure grid
  std::cout << "# create grid config" << std::endl;
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = trainDataset.getDimension();
  gridConfig.level_ = 2;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  // configure regularization
  std::cout << "# create regularization config" << std::endl;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Identity;

  //Define decomposition type
  DBMatDecompostionType dt;
  std::string decompType;
  // "LU decomposition"
  //dt = DBMatDecompLU;
  //decompType = "LU decomposition";
  // "Eigen decomposition"
  //dt = DBMatDecompEigen;
  //decompType = "Eigen decomposition";
  //"Cholesky decomposition"
  dt = DBMatDecompChol;
  decompType = "Cholesky decomposition";
  std::cout << "Decomposition type: " << decompType << std::endl;

  // if cholesky choosen -> configure adaptive refinement
  std::cout << "# create adaptive refinement config" << std::endl;
  sgpp::base::AdpativityConfiguration adaptConfig;
  adaptConfig.numRefinements_ = 2;
  adaptConfig.noPoints_ = 10;
  adaptConfig.threshold_ = 0.0;

  // initial lambda
  double lambda = 1e-1; //ToDo: set to reasonable value

  double beta = 0.0; //ToDo: set to reasonable value

  sgpp::datadriven::DBMatDensityConfiguration dconf(&gridConfig, &adaptConfig, 
                                  regularizationConfig.regType_, lambda, dt);

  // cv configuration
  double cv_lambdaStart = 1e-1;
  double cv_lambdaEnd = 1e-10;
  int cv_lambdaSteps = 5;
  bool cv_logScale = true;

  // create learner
  std::cout << "# create learner" << std::endl;
  /*sgpp::datadriven::*/LearnerSGDEOnOff learner(dconf, trainData, trainLabels, 
                                             classNum, lambda, false, 
                                             beta);

  //std::shared_ptr<sgpp::base::DataMatrix> cv_testData = std::make_shared<sgpp::base::DataMatrix>(testData);
  sgpp::base::DataMatrix* cv_testData = &testData;
  sgpp::base::DataMatrix* cv_testDataRes = nullptr; //ToDo: needed?
  learner.setCrossValidationParameters(cv_lambdaSteps, cv_lambdaStart, cv_lambdaEnd, 
                                       cv_testData, cv_testDataRes, cv_logScale);

  size_t batch_size = 1;
  unsigned int next_cv_step = 50;
  std::cout << "# start to train the learner" << std::endl;
  learner.train(batch_size, next_cv_step);

  //compute accuracy
  //ToDo:
}

/*void classifstream(ParamClass &pc) {
	std::cout << "# classifstream" << std::endl;
	pc.stest = pc.test;
	DataVector classes(0);

	clock_t end;
	clock_t begin;
	double elapsed_secs;

	std::cout << "Decomposition type: " << pc.classif_stream_decompType << std::endl;

	//Define decomposition type
	DBMatDecompostionType dt;
	if(pc.classif_stream_decompType == "LU")
	    dt = DBMatDecompLU;
	else if(pc.classif_stream_decompType == "Eigen")
	    dt = DBMatDecompEigen;
	else if(pc.classif_stream_decompType == "Chol")
	    dt = DBMatDecompChol;
	else {
	    cout << "clustc: classifstream: Unknown decomposition type " << pc.dbmat_decompType << endl;
	    exit(-1);
	}

	//Determine class lables
	vector<string> lables_fields_string;
	vector<double> lables_fields;

	tokenize(pc.classif_stream_classLables, lables_fields_string, ":");
	for(size_t i = 0; i < lables_fields_string.size(); i++){
		std::cout << "class" << i << " :" << atoi(lables_fields_string[i].c_str()) << std::endl;
	        lables_fields.push_back(atoi(lables_fields_string[i].c_str()));
	}
	
	//Set amount of classes
	int classNum = lables_fields_string.size();
	
	double lables[classNum];
	std::copy(lables_fields.begin(), lables_fields.end(), lables);

	DBMatOffline *offline;
	if(pc.classif_stream_sg_offlineFile != "") {
		std::cout << "# reading offline object..." << std::endl;
		offline = new DBMatOffline(pc.classif_stream_sg_offlineFile);
	}
	else {
		std::cout << "# creating offline object..." << std::endl;
		sgpp::base::RegularGridConfiguration gconf = {sgpp::base::GridType::Linear, pc.classif_stream_sg_dim, pc.classif_stream_sg_level };
		DBMatDensityConfiguration dconf(&gconf, NULL, sgpp::datadriven::RegularizationType::Identity, pc.classif_stream_sg_lambda, dt);
		offline = new DBMatOffline(dconf);
		offline->buildMatrix();
		begin = clock();
		offline->decomposeMatrix();
		end = clock();
		elapsed_secs = double(end-begin)/CLOCKS_PER_SEC;
		std::cout << "Decompose Mat Chol: " << elapsed_secs << std::endl;
	}

	ClassStream* prob = new ClassStream(offline, NULL, NULL, lables, classNum, pc.classif_stream_sg_lambda, false, &classes, pc.classif_stream_sg_beta);

	if(pc.classif_stream_sg_lambda < 0)
		prob->setCrossValidationParameters(pc.classif_stream_cv_lambdaStep, pc.classif_stream_cv_lambdaStart, pc.classif_stream_cv_lambdaEnd, pc.test, pc.test_cc, pc.classif_stream_cv_logscale);

	cout << "# Classification (" << prob->classname << ") " << endl;
	unsigned int next_error_step = 1;
	unsigned int next_cv_step = 0;
	if(pc.classif_stream_sg_lambda < 0)
		next_cv_step = 1;
	unsigned long long int totalInstances = 0;
	unsigned int step = 1;
	cout << "# computing L1, L2, LMax, Hellinger-Error, KL error and Memory in each statistics step." << endl;
	std::pair<DataMatrix*, DataVector*> cur_pair;
	
	std::vector<std::pair<std::list<size_t>, unsigned int> >*  RefineCoarse_ = new vector<std::pair<std::list<size_t>, unsigned int> >(prob->getNumClasses());

	//Auxiliary variables
	sgpp::base::DataVector* alpha_work;
	sgpp::base::GridIndex * gp;
	sgpp::base::DataVector p(pc.test->getNcols());

	std::list<size_t> deletedGridPoints;
	unsigned int newPoints = 0;
	
	std::vector<std::pair<DBMatOnlineDE*, double> >* onlineObjects;

	for(cur_pair = pc.input_stream->readData(pc.classif_stream_streamBlockSize); cur_pair.first != NULL; step++) {
		bool do_cv = false;
		if(next_cv_step == step) {
			do_cv = true;
			next_cv_step *= 5;
			if(step > 1) {
				//double cur_lambda = prob->getBestLambda();
				prob->setCrossValidationParameters(pc.classif_stream_cv_lambdaStep, pc.classif_stream_cv_lambdaStart, pc.classif_stream_cv_lambdaEnd, NULL, NULL, pc.classif_stream_cv_logscale);
			}
		}

		if(step == 1){
			prob->train(*(cur_pair.first), *(cur_pair.second), do_cv, RefineCoarse_);
			totalInstances += (cur_pair.first)->getNrows();
			std::cout << "Classification before applying adaptivity" << std::endl;
			classifstream_doStats(pc, prob, totalInstances, 0, step, next_error_step, false);
		}

		//Access DBMatOnlineDE-objects of all classes in order to apply adaptivity to the specific sparse grids later on
		onlineObjects = prob->getDestFunctions();

		//If the Cholesky decomposition is chosen as factorization method refinement and coarsening methods can be applied
		if(offline->getConfig()->decomp_type_== DBMatDecompChol){
			for(unsigned int i = 0; i < prob->getNumClasses(); i++){

				//Performe Coarsening
				std::cout << "Refinement and coarsening for class: " << i << std::endl;
				DBMatOnlineDE* densEst = (*onlineObjects)[i].first;
				sgpp::base::Grid* linearGitter_ = densEst->getOffline()->getGridPointer();
                                std::unique_ptr<OperationEval> opEval(sgpp::op_factory::createOperationEval(*linearGitter_));
		      		GridStorage& gridStorage = linearGitter_->getStorage();
				sgpp::base::GridGenerator& gridGen = linearGitter_->getGenerator();
		
				alpha_work = densEst->getAlpha();	
	
		      		sgpp::base::DataVector alphaWeight(alpha_work->getSize());

				//Save grid coordinates
				if(step == 1){
					classifstream_plotGrid(pc, gridStorage, step, i, deletedGridPoints, newPoints, true);
				}
				
				//Determine surpluses
		      		for (unsigned int i = 0; i < gridStorage.getSize(); i++) {
					gp = gridStorage.get(i);
					//Sets values of p to the coordinates of the given GridPoint gp
					gp->getCoords(p);
					//Multiply i-th alpha with the evaluated function at grind-point i
					alphaWeight[i] = alpha_work->get(i)*opEval->eval(*alpha_work, p);
		      		} 
				
				
				std::cout << "Size before adaptivity: " << linearGitter_->getSize(); 
				sgpp::base::HashCoarsening coarse_;
				//std::cout << std::endl << "Start coarsening" << std::endl;
		
				//Coarsening based on surpluses
				sgpp::base::SurplusCoarseningFunctor scf(alphaWeight, pc.classif_stream_number_coarse, pc.classif_stream_threshold_coarse);

				//std::cout << "Size before coarsening:" << linearGitter_->getSize() << std::endl; 
				int old_size = linearGitter_->getSize();
				coarse_.free_coarsen_NFirstOnly(linearGitter_->getStorage(), scf, alphaWeight, linearGitter_->getSize());

				//std::cout << "Size after coarsening:" << linearGitter_->getSize() << std::endl << std::endl;
				int new_size = linearGitter_->getSize();
		
				deletedGridPoints.clear();
				deletedGridPoints = coarse_.getDeletedPoints();
			
				(*RefineCoarse_)[i].first = deletedGridPoints;

				//Performe refinement
				//std::cout << "Start refinement" << std::endl;
				
				unsigned int size_before_refine = linearGitter_->getSize();

				//Refinement based on surpluses
				sgpp::base::SurplusRefinementFunctor srf(alphaWeight, pc.classif_stream_number_refine);
		      		gridGen.refine(srf);

				unsigned int size_after_refine = linearGitter_->getSize();
				std::cout << "  Size after adaptivity: " << linearGitter_->getSize() << std::endl;
				//std::cout << "Size after refinment:" << size_after_refine << std::endl << std::endl;
				newPoints = size_after_refine - size_before_refine;
				(*RefineCoarse_)[i].second = newPoints;
					
				classifstream_plotGrid(pc, gridStorage, step, i, deletedGridPoints, newPoints);

				//Apply grid changes to the Cholesky factorization	
				densEst->getOffline()->choleskyModification(newPoints, deletedGridPoints, densEst->getBestLambda());
			}
		}

		prob->train(*(cur_pair.first), *(cur_pair.second), do_cv, RefineCoarse_);

		if(step != 1){
			totalInstances += (cur_pair.first)->getNrows();
		}

		classifstream_doStats(pc, prob, totalInstances, 0, step, next_error_step, false);

		delete cur_pair.first;
		delete cur_pair.second;
		cur_pair = pc.input_stream->readData(pc.classif_stream_streamBlockSize);

	}

	step -= 1;
	classifstream_doStats(pc, prob, totalInstances, 0, step, next_error_step, true);
	std::cout << "# Stream finished after " << step << " steps and " << totalInstances << " points." << std::endl;

	delete offline;
	delete prob;
	delete RefineCoarse_;
}

inline void classifstream_predict(ParamClass &pc, ClassStream *prob, unsigned int step) {
	//prob->normalize(1000000);
	//In predict befinden sich die zugeordneten Klassen der test cases
	sgpp::base::DataVector predict = prob->predict(pc.test);
	int match, notmatch;
	match = notmatch = 0;
	unsigned int numClasses = prob->getNumClasses()+1;
	unsigned int mat[numClasses][numClasses];
	for(size_t i = 0; i < numClasses; i++)
		for(size_t j = 0; j < numClasses; j++)
			mat[i][j] = 0;
	ofstream classes;
	stringstream filename;
	filename << "/home/adrian/Dokumente/bachelorarbeit/DR10 Test/Classification/DR10_" << step << ".txt";
	classes.open(filename.str().c_str());
	for(size_t i = 0; i < pc.test->getNrows(); i++) {
		int pre = static_cast<int>(round(predict.get(i)));
		int tob = static_cast<int>(round(pc.test_cc->get(i, 0)));
		mat[pre][tob]++;
		//int what = static_cast<int>(predict.get(i))*2 + ((int)predict.get(i) == (int)pc.test_cc->get(i, 0))?1:0;
		classes << pc.test->get(i, 0) << " " << pc.test->get(i, 1) << " " << pc.test->get(i, 2) << " " << pc.test->get(i, 3) << " " << tob << " " << pre << std::endl;
		if(pre != -1 && pre == tob)
			match++;
		else
			notmatch++;
	}
	classes.close();
	for(size_t i = 0; i < numClasses; i++) {
		std::cout << std::endl << "# ";
		for(size_t j = 0; j < numClasses; j++) {
			std::cout << mat[i][j] << " ";
		}
	}
	std::cout << std::endl << " Match: " << match << ", not match: " << notmatch << ", rate: " << (double)(((double)match*100.)/(double)(match+notmatch));
	return;
}

inline void classifstream_doStats(ParamClass &pc, ClassStream *prob, unsigned long long int totalInstances, int start_mem_value, unsigned int step, unsigned int &next_error_step, bool last) {
	if(step == next_error_step || last) {
		cout << totalInstances << " " << step << " ";
		if(pc.classif_stream_sg_lambda < 0)
			cout << prob->getBestLambda();
		else
			cout << pc.classif_stream_sg_lambda;
		//int now_mem_value = getMemValue();
		//std::cout << " Memory:" << now_mem_value << " ";
		next_error_step *= 2;
		//destream_doError(pc, prob);
		classifstream_predict(pc, prob, step);
		std::cout << std::endl;
	}
}

inline void tokenize(const string& str, vector<string>& tokens, const string& delimiters) {
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

inline int getMemValue(){ //Note: this value is in KB!
	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];
	
	
	while (fgets(line, 128, file) != NULL){
		if (strncmp(line, "VmSize:", 7) == 0){
			result = parsekLine(line);
			break;
		}
	}
	fclose(file);
	return result;
}

inline int parsekLine(char* line){
	int i = strlen(line);
	while (*line < '0' || *line > '9') line++;
	line[i-3] = '\0';
	i = atoi(line);
	return i;
}

inline void classifstream_plotGrid(ParamClass &pc, GridStorage& gridStorage, size_t step, size_t klasse, std::list<size_t> deletedPoints, int newPoints, bool first){

	sgpp::base::GridIndex * gp;
	string fileName;	
	ofstream grids;
	string string1 = "/home/adrian/Dokumente/bachelorarbeit/DR10 Test/Gitter/DR10_Gitter";
	string string2 = ".txt";
	fileName = string1 + std::to_string(step) + std::to_string(klasse) + string2;
	if(first){
		string stringFirst = "first";
		fileName = string1 + std::to_string(step) + std::to_string(klasse) + stringFirst + string2;
	}
	grids.open(fileName.c_str());
	sgpp::base::DataVector p(pc.test->getNcols());

	int counter = 1;
	if(deletedPoints.size() != 0){
		for (std::list<size_t>::iterator it = deletedPoints.begin(); it != deletedPoints.end(); it++){
			if(counter < deletedPoints.size())
				grids << *it << " ";
			else
				grids << *it << " " << newPoints << std::endl;
			counter++;
		}
	}else{
		grids << newPoints << std::endl;
	}
	for (unsigned int i = 0; i < gridStorage.getSize(); i++) {
		gp = gridStorage.get(i);
		//Sets values of p to the coordinates of the given GridPoint gp
		gp->getCoords(p);
		for(size_t j = 0; j < pc.classif_stream_sg_dim; j++) {
			if(j < (pc.classif_stream_sg_dim - 1)){
				grids << p[j] << " ";
			}else{
				grids << p[j] << std::endl;
			}
		}
	} 
	grids.close();
	return;

}

//Compare two DataMatrix-objects
bool compareDataMatrix(DataMatrix* mat1, DataMatrix* mat2){
	
	//std::cout << "Compare rein" <<std::endl;
	unsigned int size_1 = mat1->getSize();
	unsigned int size_2 = mat2->getSize();
	
	if((mat1->getNcols() != mat2->getNcols()) || (mat1->getNrows() != mat2->getNrows())){
		std::cout << "Both matrices need to have same size" << std::endl;
		return false;
	}
	
	sgpp::base::DataMatrix* temp_copy = new sgpp::base::DataMatrix(mat1->getNrows(), mat1->getNcols(), 0.0);
	temp_copy->copyFrom(*mat1);
	temp_copy->sub(*mat2);
	
	temp_copy->abs();
	float max = temp_copy->max();

	delete temp_copy;

	if( fabs(max) > std::numeric_limits<float>::epsilon()){
		return false;
	} else {
		return true;
	}
	
	return 0;
}*/



