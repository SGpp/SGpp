// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSChol.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
//#include <sgpp/base/grid/type/LinearTrapezoidBoundaryGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
//#include <sgpp/base/operation/OperationMatrix.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
//#include <sgpp/combigrid/operation/hash/OperationMatrixLTwoExplicitFullGrid.hpp>

#ifdef USE_GSL
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#endif /* USE_GSL */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <ctime>

using namespace std;

DBMatOffline::DBMatOffline(sgpp::datadriven::DBMatDensityConfiguration& oc) :
		config_(&oc), lhsMatrix_(NULL), constructed_(false), decomposed_(false),
                ownsConfig_(false), perm_(NULL), grid_(NULL) {
}

DBMatOffline::DBMatOffline () :
		config_(NULL), lhsMatrix_(NULL), constructed_(false), decomposed_(false), ownsConfig_(
						false), perm_(NULL), grid_(NULL) {
		}

DBMatOffline::DBMatOffline(string fname) :
		config_(NULL), lhsMatrix_(NULL), constructed_(false), decomposed_(
				false), ownsConfig_(true), perm_(NULL), grid_(NULL) {

	//std::cout << "START READING MATRIX" << std::endl;
	//sgpp::base::SGppStopwatch* myStopwatch = new sgpp::base::SGppStopwatch();
	//myStopwatch->start();

	//Read configuration
	FILE * f = fopen(fname.c_str(), "rb");
	string line("");
	if (f == NULL) {
		cout << "DBMatOffline: Error opening file " << line << endl;
		exit(-1);
	}
	char c = fgetc(f);
	line += c;
	while (c != '\n') {
		c = fgetc(f);
		line += c;
	}

	vector<string> tokens;
	string delim(",");
	Tokenize(line, tokens, delim);

        // ToDo: full grid not supported
	bool fullgrid = atoi(tokens[0].c_str());
	if ((fullgrid && tokens.size() != 8)
			|| (!fullgrid && tokens.size() != 7)) {
		cout << "DBMatOffline: Wrong file format: " << fname.c_str()
				<< endl;
		exit(-1);
	}

	sgpp::base::GridType grid_type = (sgpp::base::GridType) atoi(tokens[1].c_str());
	
	int grid_dim = atoi(tokens[2].c_str());
	int grid_level = atoi(tokens[3].c_str());
	sgpp::datadriven::RegularizationType reg =
			(sgpp::datadriven::RegularizationType) atoi(tokens[4].c_str());
	double lambda = atof(tokens[5].c_str());
	DBMatDecompostionType decomp = (DBMatDecompostionType) atoi(
			tokens[6].c_str());

	sgpp::base::RegularGridConfiguration gconf = { grid_type, grid_dim,
			grid_level };
	config_ = new sgpp::datadriven::DBMatDensityConfiguration(&gconf, NULL, reg, lambda,
			decomp);

	int size;

	//Build grid
	InitializeGrid();

	// check if grid was created
	if (grid_ == NULL)
		return;

	size = grid_->getStorage().getSize(); //Size of the (quadratic) matrices A and C
	
	//Read matrix
#ifdef USE_GSL
	gsl_matrix* a;
	if (decomp == DBMatDecompLU) {
		a = gsl_matrix_alloc(size, size);
	} else if (decomp == DBMatDecompEigen) {
		a = gsl_matrix_alloc(size + 1, size);
	} else if (decomp == DBMatDecompChol){
		a = gsl_matrix_alloc(size, size);
	} else {
		throw sgpp::base::application_exception(
				"Unsupported decomposition type!");
	}

	gsl_matrix_fread(f, a);
	lhsMatrix_ = new sgpp::base::DataMatrix(a->data, a->size1, a->size2);
	//Read permutation
	if (decomp == DBMatDecompLU) {
		perm_ = gsl_permutation_alloc(size);
		gsl_permutation_fread(f, perm_);
	}

	fclose(f);
	gsl_matrix_free(a);

	constructed_ = true;
	decomposed_ = true;

	//std::cout << "Time: " << myStopwatch->stop() << std::endl;
#endif /* USE_GSL */
}

DBMatOffline::DBMatOffline(const DBMatOffline &old) :
		config_(NULL), lhsMatrix_(NULL), constructed_(true), decomposed_(
				true), ownsConfig_(true), perm_(NULL), grid_(NULL) {

	if (!old.decomposed_) {
		throw sgpp::base::application_exception("Matrix not yet decomposed");
		return;
	}

	//Copy configuration	
	sgpp::base::RegularGridConfiguration gconf = {old.config_->grid_type_, old.config_->grid_dim_,
			                               old.config_->grid_level_};

        sgpp::base::AdpativityConfiguration adaptConfig;
        adaptConfig.numRefinements_ = old.config_->numRefinements_;
        adaptConfig.threshold_ = old.config_->ref_threshold_;
        adaptConfig.noPoints_ = old.config_->ref_noPoints_;

	config_ = new sgpp::datadriven::DBMatDensityConfiguration(&gconf, &adaptConfig, old.config_->regularization_, old.config_->lambda_,
			old.config_->decomp_type_);
	
	//Copy Matrix
	lhsMatrix_ = new sgpp::base::DataMatrix(*(old.lhsMatrix_));
#ifdef USE_GSL
	//Copy Permutation (if existing)
	if (config_->decomp_type_ == DBMatDecompLU) {
		perm_ = gsl_permutation_alloc(old.grid_->getStorage().getSize());
		gsl_permutation_memcpy(perm_, old.perm_);
	}
#endif /* USE_GSL */
	//Initialize Grid
	InitializeGrid();
}

DBMatOffline::~DBMatOffline() {
	if (lhsMatrix_ != NULL)
		delete lhsMatrix_;
	if (perm_ != NULL)
		delete perm_;
	if (grid_ != NULL)
		delete grid_;
	/*if (fullgrid_ != NULL)
		delete fullgrid_;*/
	if (ownsConfig_ && config_ != NULL)
		delete config_;
}

sgpp::datadriven::DBMatDensityConfiguration* DBMatOffline::getConfig() {
	return config_;
}

sgpp::base::DataMatrix * DBMatOffline::getDecomposedMatrix() {
	//if (decomposed_) {
		return lhsMatrix_;
	/*} else {
		throw sgpp::base::data_exception("Matrix was not decomposed, yet!");
	}*/
}

//combigrid::FullGrid<double>& DBMatOffline::getFullGrid() {
//	return *fullgrid_;
//}

sgpp::base::Grid & DBMatOffline::getGrid() {
	return *grid_;
}

sgpp::base::Grid* DBMatOffline::getGridPointer() {
  return grid_;
}

void DBMatOffline::permuteVector(sgpp::base::DataVector& b) {
#ifdef USE_GSL
	if (decomposed_) {
		if (config_->decomp_type_ == DBMatDecompLU) {
			gsl_permute(perm_->data, b.getPointer(), 1, b.getSize());
		} else {
			//No permutation needed for other decomposition types!
		}
	} else {
		throw sgpp::base::data_exception("Matrix was not decomposed, yet!");
	}
#endif /* USE_GSL */
}

void DBMatOffline::InitializeGrid() {
	//if (config_->grid_type_ == sgpp::base::GridType::LinearTrapezoidBoundary) {
	//	grid_ = new sgpp::base::LinearTrapezoidBoundaryGrid(config_->grid_dim_);
	//} else
	if (config_->grid_type_ == sgpp::base::GridType::ModLinear) {
		grid_ = new sgpp::base::ModLinearGrid(config_->grid_dim_);
	} else if (config_->grid_type_ == sgpp::base::GridType::Linear) {
		grid_ = new sgpp::base::LinearGrid(config_->grid_dim_);
	} else {
		grid_ = NULL;
		throw sgpp::base::application_exception(
				"LearnerBase::InitializeGrid: An unsupported grid type was chosen!");
	}

// Generate regular Grid with LEVELS Levels
	grid_->getGenerator().regular(config_->grid_level_);
}

void DBMatOffline::buildMatrix() {
	if (constructed_) { //Already constructed, do nothing
		return;
	}

	unsigned int size;

	InitializeGrid();

        // check if grid was created
	if (grid_ == NULL)
		return;

	size = grid_->getStorage().getSize(); //Size of the (quadratic) matrices A and C
	


//Construct matrix A
	lhsMatrix_ = new sgpp::base::DataMatrix(size, size);
	
	std::unique_ptr<sgpp::base::OperationMatrix> op(sgpp::op_factory::createOperationLTwoDotExplicit(lhsMatrix_, *grid_));
	
	if (config_->decomp_type_ == DBMatDecompLU || config_->decomp_type_ == DBMatDecompChol) { //Add regularization term:

//Construct matrix lambda * C (just use identity for C)
		sgpp::base::DataMatrix lambdaC(size, size);
		if (config_->regularization_ == sgpp::datadriven::RegularizationType::Identity) {
			lambdaC.setAll(0.);
			for (unsigned int i = 0; i < size; i++) {
				lambdaC.set(i, i, config_->lambda_);
			}
		} else {
			throw sgpp::base::operation_exception(
					"Unsupported regularization type");
		}

//Compute A + lambda * C:

		lhsMatrix_->add(lambdaC);
	}

	constructed_ = true;
}

void DBMatOffline::decomposeMatrix() {
#ifdef USE_GSL
	if (constructed_) {
		if (decomposed_) //Already decomposed => Do nothing
		{
			return;
		} else {
			unsigned int n = lhsMatrix_->getNrows();
			if (config_->decomp_type_ == DBMatDecompLU) {
				gsl_matrix_view m = gsl_matrix_view_array(
						lhsMatrix_->getPointer(), n, n); //Create GSL matrix view for decomposition
				int* sig = (int*) malloc(sizeof(int));
				perm_ = gsl_permutation_alloc(n); //allocate permutation

				gsl_linalg_LU_decomp(&m.matrix, perm_, sig);

				delete sig;

			} else if (config_->decomp_type_ == DBMatDecompEigen) {
				gsl_matrix_view m = gsl_matrix_view_array(
						lhsMatrix_->getPointer(), n, n); //Create GSL matrix view for decomposition

				gsl_matrix* q = gsl_matrix_alloc(n, n); //Stores the eigenvectors
				gsl_vector* e = gsl_vector_alloc(n); //Stores the eigenvalues

				gsl_eigen_symmv_workspace* ws = gsl_eigen_symmv_alloc(n);
				gsl_eigen_symmv(&m.matrix, e, q, ws);
				gsl_eigen_symmv_free(ws);

				//Create an (n+1)*n matrix to store eigenvalues and -vectors:
				delete lhsMatrix_;
				lhsMatrix_ = new sgpp::base::DataMatrix(n + 1, n);

				for (unsigned int r = 0; r < n; r++) {
					for (unsigned int c = 0; c < n; c++) {
						lhsMatrix_->set(r, c, gsl_matrix_get(q, r, c));
					}
				}
				for (unsigned int c = 0; c < n; c++) {
					lhsMatrix_->set(n, c, gsl_vector_get(e, c));
				}

				delete e;
				delete q;
			} else if (config_->decomp_type_ == DBMatDecompChol){

				gsl_matrix_view m = gsl_matrix_view_array(
							lhsMatrix_->getPointer(), n, n); //Create GSL matrix view for decomposition
				//Perform Cholesky decomposition
				gsl_linalg_cholesky_decomp(&m.matrix);

				//Isolate lower triangular matrix
				for(unsigned int i = 0; i < n; i++){
					for(unsigned int j = 0; j < n; j++){
						if(i < j){
						   lhsMatrix_->set(i,j,0);
						}
					}
				}
			} else {
				//No other decomposition types, yet
				return;
			}
			decomposed_ = true;
		}
	} else {
		throw sgpp::base::data_exception(
				"Matrix has to be constructed before it can be decomposed!");
	}
#endif /* USE_GSL */
}

void DBMatOffline::store(string fileName) {
	if (!decomposed_) {
		throw sgpp::base::application_exception("Matrix not yet decomposed");
		return;
	}

	//Write configuration
	FILE * f = fopen(fileName.c_str(), "w");
	if (f <= 0) {
		cout << "libtool: DBMatOffline: Error opening file " << fileName.c_str()
				<< endl;
		exit(-1);
	}

	fprintf(f, "%d,%d,%d,%d,%d,%.10e,%d", config_->grid_type_,
			config_->grid_dim_, config_->grid_level_, config_->regularization_,
			config_->lambda_, config_->decomp_type_);
	fprintf(f, "\n");
	fclose(f);
#ifdef USE_GSL
	//Write Matrix
	f = fopen(fileName.c_str(), "ab");
	gsl_matrix_view m = gsl_matrix_view_array(lhsMatrix_->getPointer(),
			lhsMatrix_->getNrows(), lhsMatrix_->getNcols());
	gsl_matrix_fwrite(f, &m.matrix);

	//Write Permutation (if existing)
	if (config_->decomp_type_ == DBMatDecompLU) {
		gsl_permutation_fwrite(f, perm_);
	}

	fclose(f);
#endif /* USE_GSL */
}

void DBMatOffline::printMatrix() {
	std::cout << "Size: " <<  lhsMatrix_->getNrows() << " , " << lhsMatrix_->getNcols() << std::endl;
	for (unsigned int r = 0; r < lhsMatrix_->getNrows(); r++) {
		for (unsigned int c = 0; c < lhsMatrix_->getNcols(); c++) {
			std::cout << lhsMatrix_->get(r, c) << ", ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


void DBMatOffline::Tokenize(string& str, vector<string>& tokens,
		string& delimiters) {
	if (!strcmp(delimiters.c_str(), "")) {
		delimiters = " ";
	}
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos) {
	// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
	// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
	// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}

}

void DBMatOffline::choleskyModification(int newPoints,std::list<size_t> deletedPoints, double lambda){
	//Start coarsening
		
	//If list 'deletedPoints' is not empty, grid points got removed
	if(deletedPoints.size() > 0){

		int new_size = grid_->getSize();
		int old_size = new_size - newPoints + deletedPoints.size();
	
		//'c' is the threshold to decide what kind of Cholesky modification should be performed to take
		//account for the removed grid points
		double c = old_size/10;
		//Indicates which permutation to perform
		bool start1 = false;
		size_t index_coarse = 0;

		//If index of removed grid point <= c-> Permutation towards first row/column -> job = 1
		size_t coarseCount_1 = 0;
		//If index of removed grid point > c -> Permutation towards last row/column -> job = 2
		size_t coarseCount_2 = 0;
		
	  	for (std::list<size_t>::reverse_iterator it = deletedPoints.rbegin(); it != deletedPoints.rend(); it++){
			
			//Indicates current row/column to permute
			index_coarse = *it + 1 + coarseCount_1;

			if(index_coarse > c && start1 == false){
			//Is accessed if index is larger than 'c'
				choleskyPermutation(index_coarse, old_size - coarseCount_2, 2);
				coarseCount_2++;
			}else if(start1 == false){
			//Is accessed the first 'index_coarse' is less than 'c'
				start1 = true;
				//Delete last 'coarseCount_2' rows/columns of current Cholesky factor
				lhsMatrix_->resizeQuadratic(old_size - coarseCount_2);
			}

			if(start1 == true){
			//Is accessed if 'index_coarse' is less than 'c'
				choleskyPermutation(1 + coarseCount_1, index_coarse, 1);
				coarseCount_1++;
			}
		}

		if(coarseCount_1 > 0){
		//If some indices have been less than 'c'
			sgpp::base::DataMatrix* update_matrix = new sgpp::base::DataMatrix(lhsMatrix_->getNrows(), lhsMatrix_->getNcols());
			update_matrix->copyFrom(*lhsMatrix_);
			//Resize copy of current Cholesky factor to receive a matrix of rank one update vectors
			update_matrix->resizeToSubMatrix(coarseCount_1 + 1, 1, lhsMatrix_->getNrows(), coarseCount_1);
			//Resize current Cholesky factor to required submatrix for necessary rank one updates
			lhsMatrix_->resizeToSubMatrix(coarseCount_1 + 1, coarseCount_1 + 1, lhsMatrix_->getNrows(), lhsMatrix_->getNrows());
			sgpp::base::DataVector* temp_col = new sgpp::base::DataVector(update_matrix->getNrows());

			//'coarseCount_1' many rank one updates based on the columns of 'update_matrix' are performed
			DBMatDMSChol cholsolver;
			for(size_t i = 0; i < coarseCount_1; i++){
				update_matrix->getColumn(i, *temp_col);
				cholsolver.choleskyUpdate(*lhsMatrix_, temp_col, false);
			}
			delete update_matrix;
			delete temp_col;
		}else{
		//If no indices have been less than 'c'
			lhsMatrix_->resizeQuadratic(old_size - coarseCount_2);
		}
	}

	//Start refinement
	if(newPoints > 0){
		
		size_t gridSize = grid_->getStorage().getSize();
	      	size_t gridDim  = grid_->getStorage().getDimension();
		
		//DataMatrix to collect vectors to append 
		sgpp::base::DataMatrix* mat_refine = new sgpp::base::DataMatrix(gridSize, newPoints);

	      	sgpp::base::DataMatrix level(gridSize, gridDim);
	      	sgpp::base::DataMatrix index(gridSize, gridDim);
		
		grid_->getStorage().getLevelIndexArraysForEval(level, index);
		float lambda_conf = lambda;
		//Loop to calculate all L2-products of added points based on the hat-function as basis function
		for (size_t i = 0; i < gridSize; i++) {
			for (size_t j = gridSize-newPoints; j < gridSize; j++) {
				  float res = 1;

				  for (size_t k = 0; k < gridDim; k++) {
				    float lik = level.get(i, k);
				    float ljk = level.get(j, k);
				    float iik = index.get(i, k);
				    float ijk = index.get(j, k);

				    if (lik == ljk) {
				      if (iik == ijk) {
					//Use formula for identical ansatz functions:
					res *= 2 / lik / 3;
				      } else {
					//Different index, but same level => ansatz functions do not overlap:
					res = 0.;
					break;
				      }
				    } else {
				      if (std::max((iik - 1) / lik, (ijk - 1) / ljk)
					  >= std::min((iik + 1) / lik, (ijk + 1) / ljk)) {
					//Ansatz functions do not not overlap:
					res = 0.;
					break;
				      } else {
					//Use formula for different overlapping ansatz functions:
					if (lik > ljk) { //Phi_i_k is the "smaller" ansatz function
					  float diff = (iik / lik) - (ijk / ljk); // x_i_k - x_j_k
					  float temp_res = fabs(diff - (1 / lik))
						             + fabs(diff + (1 / lik)) - fabs(diff);
					  temp_res *= ljk;
					  temp_res = (1 - temp_res) / lik;
					  res *= temp_res;
					} else { //Phi_j_k is the "smaller" ansatz function
					  float diff = (ijk / ljk) - (iik / lik); // x_j_k - x_i_k
					  float temp_res = fabs(diff - (1 / ljk))
						             + fabs(diff + (1 / ljk)) - fabs(diff);
					  temp_res *= lik;
					  temp_res = (1 - temp_res) / ljk;
					  res *= temp_res;
					}
				      }
				    }
				  }
				  //The new Rows/Cols are stored in mat_refine

				  //add current lambda to lower diagonal elements of mat_refine
				  if(i == j){
					mat_refine->set(i, j - gridSize + newPoints, res + lambda_conf);
				  }else{
				  	mat_refine->set(i, j - gridSize + newPoints , res);
				  }
			}
		}
			
		//Resize Cholesky factor to new 'gridSize' before 'choleskyAddPoint' is applied
		//in order to save runtime
		this->lhsMatrix_->resizeQuadratic(gridSize);

		//Now its time to call 'choleskyAddPoint''countNewGridPoints' often
		sgpp::base::DataVector* temp_col = new sgpp::base::DataVector(gridSize);
		for(int j = gridSize-newPoints; j < gridSize; j++){
			temp_col->resizeZero(gridSize);
			mat_refine->getColumn(j - gridSize + newPoints, *temp_col);	
			temp_col->resizeZero(j + 1);
			choleskyAddPoint(temp_col, j);
		}

		delete temp_col;
	}

	return;
}

void DBMatOffline::choleskyAddPoint(sgpp::base::DataVector* newCol, size_t size){
#ifdef USE_GSL	
	if (!decomposed_) {
		throw sgpp::base::data_exception("Matrix was not decomposed, yet!");
	}

	sgpp::base::DataMatrix* mat_ = lhsMatrix_;
	//Size of provided memory for Cholesky factor, because the allocations take place in 'choleskyModifications'
	unsigned int size_full = mat_->getNrows();
	//Size of Cholesky factor after adding 'newCol'
	unsigned int size_up = newCol->getSize();
	
	if(size_up != (size+1)){
		throw sgpp::base::data_exception("Size of update vector newCol needs to be 1 dim larger than the underlying decomposed matrix!");
	}
	
	//Create GSL matrix view for update procedures
	gsl_matrix_view m_full = gsl_matrix_view_array(mat_->getPointer(), size_full, size_full); 
	//Access submatrx since mat_ has already expanded to save alocations procedures
	gsl_matrix_view m = gsl_matrix_submatrix(&m_full.matrix, 0, 0, size, size);

	gsl_vector_view vvec = gsl_vector_view_array(newCol->getPointer(), size_up);
	gsl_vector*  wkvec_a = gsl_vector_calloc(size);

	//Extract newCol(0:size-1) = a
	gsl_vector_view c = gsl_vector_subvector(&vvec.vector, 0, size);

	//Solve system a = Lc
	gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, &m.matrix, &c.vector);
	gsl_blas_dcopy(&c.vector, wkvec_a);

	double phi;
	gsl_blas_ddot (wkvec_a, &c.vector, &phi);
	
	//Compute d = sqrt(newCol[last] - phi)
	double last = gsl_vector_get(&vvec.vector, size_up - 1);
	phi = last - phi;
	//Ensure 'phi' is larger 0
	if(phi <= 0){
		throw sgpp::base::data_exception("Resulting matrix is at least not numerical positive definite!");
	}else{
		phi = sqrt(phi);
	}
	
	//Modify Choleskyfactor m (L)

	//DataVector full of zeros
	sgpp::base::DataVector* zeros = new sgpp::base::DataVector(size_full, 0.0);

	//Add 'newCol' and 'zeros' to current Cholesky factor 'lhsMatrix_'
	mat_->setColumn(size, *zeros);
	newCol->set(size_up-1, phi);
	newCol->resizeZero(size_full);
	mat_->setRow(size, *newCol);

	gsl_vector_free(wkvec_a);
	delete zeros;

	return;
#endif /* USE_GSL */
}

void DBMatOffline::choleskyPermutation(int k, int l, int job){
#ifdef USE_GSL	
	if (!decomposed_) {
		throw sgpp::base::data_exception("Matrix was not decomposed, yet!");
	}
	
	sgpp::base::DataMatrix* mat_ = lhsMatrix_;
	int size = mat_->getNrows();

	if(k > l){
		throw sgpp::base::data_exception("l needs to be larger than k");
	}else if(l > size){
		throw sgpp::base::data_exception("l needs to be smaller than the column size of the matrix");
	}else if(l == k){
		return;
	}
	
	//Determine amount of necessary Givens rotations and rows to swap
	int count_perm = l-k;
	
	gsl_matrix_view m = gsl_matrix_view_array(mat_->getPointer(), size, size); //Create GSL matrix view for update procedures
		

	//Define and declare Workingvector, Cosine- and Sinevector
	gsl_vector*  svec  = gsl_vector_calloc(count_perm);
	gsl_vector*  cvec  = gsl_vector_calloc(count_perm);
	
	double* givens_zero;
	double* tbuff = m.matrix.data;
	
	if(job == 2){
		//Permute upper triangular - job = 2        => left circular shift 
		//1,...,k-1,k,k+1, ..., l-1,l,l+1, ..,size  => 1,...,k-1,k+1, ..., l-1,l,k,l+1,..., size
		for(int i = k - 1; i < l - 1; i++){
			gsl_matrix_swap_rows(&m.matrix, i, i + 1);
		}

		tbuff += (k-1)*size + (k-1);
		for(int j = 1; j <= count_perm; j++){

			givens_zero = gsl_matrix_ptr (&m.matrix, k + j - 2, k + j - 1);
			//Givensrotation in (k+i-1,k+i)-Plane
			gsl_blas_drotg(tbuff, givens_zero, cvec->data+j-1, svec->data+j-1);
			//Access columns to modify via Givens rotation
			gsl_vector_view diag_sub     = gsl_matrix_subcolumn(&m.matrix, k + j - 2, k + j - 1, size - k - j + 1);
			gsl_vector_view low_diag_sub = gsl_matrix_subcolumn(&m.matrix, k + j - 1, k + j - 1, size - k - j + 1);
			gsl_matrix_set (&m.matrix, k + j - 2, k + j - 1, 0.0);
			//Apply Givens rotation	
			gsl_blas_drot(&diag_sub.vector, &low_diag_sub.vector, cvec->data[j-1], svec->data[j-1]);
			tbuff += (size + 1);
			
		}	
	}else if(job == 1){
		//Permute upper triangular - job = 1       => right circular shift 
		//1,...,k-1,k,k+1, ..., l-1,l,l+1,...size  => 1,...,k-1,l,k,k+1, ..., l-1,l+1,...size 
		for(int i = l - 1; i > k - 1; i--){
			gsl_matrix_swap_rows(&m.matrix, i, i - 1);
		}
		
		tbuff += (k-1) * size + (l -2);
		for(int j = 1; j <= count_perm; j++){
			
			givens_zero = gsl_matrix_ptr (&m.matrix, k - 1, l - j);
			//Givensrotation in (l-i,l-i+1)-Plane
			gsl_blas_drotg(tbuff, givens_zero, cvec->data+j-1, svec->data+j-1);
			//Access columns to modify via Givens rotation
			gsl_vector_view diag_sub     = gsl_matrix_subcolumn(&m.matrix, l - j - 1, l - j, size - l + j);
			gsl_vector_view low_diag_sub = gsl_matrix_subcolumn(&m.matrix, l - j, l - j , size - l + j);
			gsl_matrix_set (&m.matrix, k - 1, l - j, 0.0);
			//Apply Givens rotation	
			gsl_blas_drot(&diag_sub.vector, &low_diag_sub.vector, cvec->data[j-1], svec->data[j-1]);
			tbuff -= 1;
		}
	}
	
	gsl_vector_free(svec);
	gsl_vector_free(cvec);

	return;	
#endif /* USE_GSL */
}
