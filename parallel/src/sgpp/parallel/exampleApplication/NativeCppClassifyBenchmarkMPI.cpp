// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <mpi.h>

#include <omp.h>

#include <sgpp_mpi.hpp>
#include <sgpp_base.hpp>
#include <sgpp_parallel.hpp>
#include <sgpp_datadriven.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/tools/DatasetGenerator.hpp>

#include <string>
#include <iostream>
#include <ostream>
#include <cstdlib>
#include <fstream>

// print grid in gnuplot readable format (1D and 2D only)
//#define GNUPLOT
//#define GRDIRESOLUTION 100

bool bUseFloat;
bool bUseRecursion;

std::string ggridtype;
std::string gdataFile;
std::string gtestFile;
bool gisRegression;
SGPP::datadriven::ClassificatorQuality gTrainQual;
SGPP::datadriven::ClassificatorQuality gTestQual;
SGPP::datadriven::LearnerTiming gtimings;
double gtrainAcc;
double gtestAcc;
SGPP::solver::SLESolverConfiguration gSLEfinal;
SGPP::base::AdpativityConfiguration gAdapConfig;
int gstart_level;
double glambda;

//void storeROCcurve(SGPP::base::DataMatrix& ROC_curve, std::string tFilename)
//{
//  std::ofstream fileout;
//
//  // Open filehandle
//  fileout.open(tFilename.c_str());
//
//  // plot values
//  for (size_t i = 0; i < ROC_curve.getNrows(); i++)
//  {
//    fileout <<  ROC_curve.get(i, 0) << " " << ROC_curve.get(i, 1) << std::endl;
//  }
//
//  // close filehandle
//  fileout.close();
//}

void printSettings(std::string dataFile, std::string testFile, bool isRegression,
                   const SGPP::base::RegularGridConfiguration& GridConfig, const SGPP::solver::SLESolverConfiguration& SolverConfigRefine,
                   const SGPP::solver::SLESolverConfiguration& SolverConfigFinal, const SGPP::base::AdpativityConfiguration& AdaptConfig,
                   const double lambda, const SGPP::parallel::VectorizationType vecType) {
  if (SGPP::parallel::myGlobalMPIComm->getMyRank() != 0) {
    return;
  }

  std::cout << std::endl;

  std::cout << "Train dataset: " << dataFile << std::endl;
  std::cout << "Test dataset: " << testFile << std::endl;
  std::cout << "Startlevel: " << GridConfig.level_ << std::endl << std::endl;

  std::cout << "Num. Refinements: " << AdaptConfig.numRefinements_ << std::endl;
  std::cout << "Refine Threshold: " << AdaptConfig.threshold_ << std::endl;
  std::cout << "Refine number points: " << AdaptConfig.noPoints_ << std::endl << std::endl;

  std::cout << "Max. CG Iterations (refine): " << SolverConfigRefine.maxIterations_ << std::endl;
  std::cout << "CG epsilon (refine): " << SolverConfigRefine.eps_ << std::endl;
  std::cout << "Max. CG Iterations (final): " << SolverConfigFinal.maxIterations_ << std::endl;
  std::cout << "CG epsilon (final): " << SolverConfigFinal.eps_ << std::endl << std::endl;

  std::cout << "Lambda: " << lambda << std::endl << std::endl;

  if (bUseFloat) {
    std::cout << "Precision: Single Precision (float)" << std::endl << std::endl;
  } else {
    std::cout << "Precision: Double Precision (double)" << std::endl << std::endl;
  }

  if (vecType == SGPP::parallel::X86SIMD) {
#if defined(__SSE3__) && !defined(__AVX__)
    std::cout << "Vectorized: X86SIMD (SSE3)" << std::endl << std::endl;
#endif
#if defined(__SSE3__) && defined(__AVX__)
    std::cout << "Vectorized: X86SIMD (AVX)" << std::endl << std::endl;
#endif
  } else if (vecType == SGPP::parallel::OpenCL) {
    std::cout << "Vectorized: OpenCL (NVIDIA Fermi optimized)" << std::endl << std::endl;
  } else if (vecType == SGPP::parallel::Hybrid_X86SIMD_OpenCL) {
#if defined(__SSE3__) && !defined(__AVX__)
    std::cout << "Vectorized: Hybrid, SSE3 and OpenCL (NVIDIA Fermi optimized)" << std::endl << std::endl;
#endif
#if defined(__SSE3__) && defined(__AVX__)
    std::cout << "Vectorized: Hybrid, AVX and OpenCL (NVIDIA Fermi optimized)" << std::endl << std::endl;
#endif
  } else if (vecType == SGPP::parallel::ArBB) {
    std::cout << "Vectorized: Intel Array Building Blocks" << std::endl << std::endl;
  } else if (vecType == SGPP::parallel::MIC) {
    std::cout << "Vectorized: Intel MIC Architecture" << std::endl << std::endl;
  } else if (vecType == SGPP::parallel::Hybrid_X86SIMD_MIC) {
#if defined(__SSE3__) && !defined(__AVX__)
    std::cout << "Vectorized: Hybrid, SSE3 and Intel MIC Architecture" << std::endl << std::endl;
#endif
#if defined(__SSE3__) && defined(__AVX__)
    std::cout << "Vectorized: Hybrid, AVX and Intel MIC Architecture" << std::endl << std::endl;
#endif
  } else {
    std::cout << "Scalar Version" << std::endl << std::endl;
  }

  if (isRegression) {
    std::cout << "Mode: Regression" << std::endl << std::endl;
  } else {
    std::cout << "Mode: Classification" << std::endl << std::endl;
  }

  if (GridConfig.type_ == SGPP::base::Linear) {
    std::cout << "chosen gridtype: Linear" << std::endl << std::endl;
  } else if (GridConfig.type_ == SGPP::base::LinearBoundary) {
    std::cout << "chosen gridtype: LinearBoundary" << std::endl << std::endl;
  } else {
    const char* modlinear_mode = getenv("SGPP_MODLINEAR_EVAL");

    if (modlinear_mode == NULL) {
      modlinear_mode = "mask";
    }

    std::cout << "chosen gridtype: ModLinear (" << modlinear_mode << ")" << std::endl << std::endl;
  }


}

void printResults() {
  if (SGPP::parallel::myGlobalMPIComm->getMyRank() != 0) {
    return;
  }

  if (gisRegression) {
    std::cout << "training MSE: " << gtrainAcc << std::endl;
    std::cout << "testing MSE: " << gtestAcc << std::endl;
  } else {
    double trainSens = static_cast<double>(gTrainQual.truePositive_) / static_cast<double>(gTrainQual.truePositive_ + gTrainQual.falseNegative_);
    double trainSpec = static_cast<double>(gTrainQual.trueNegative_) / static_cast<double>(gTrainQual.trueNegative_ + gTrainQual.falsePositive_);
    double trainPrec = static_cast<double>(gTrainQual.truePositive_) / static_cast<double>(gTrainQual.truePositive_ + gTrainQual.falsePositive_);

    double testSens = static_cast<double>(gTestQual.truePositive_) / static_cast<double>(gTestQual.truePositive_ + gTestQual.falseNegative_);
    double testSpec = static_cast<double>(gTestQual.trueNegative_) / static_cast<double>(gTestQual.trueNegative_ + gTestQual.falsePositive_);
    double testPrec = static_cast<double>(gTestQual.truePositive_) / static_cast<double>(gTestQual.truePositive_ + gTestQual.falsePositive_);

    std::cout << "training accuracy: " << gtrainAcc << std::endl;
    std::cout << "training sensitivity: " << trainSens << std::endl;
    std::cout << "training specificity: " << trainSpec << std::endl;
    std::cout << "training precision: " << trainPrec << std::endl << std::endl;
    std::cout << "training true positives: " << gTrainQual.truePositive_ << std::endl;
    std::cout << "training true negatives: " << gTrainQual.trueNegative_ << std::endl;
    std::cout << "training false positives: " << gTrainQual.falsePositive_ << std::endl;
    std::cout << "training false negatives: " << gTrainQual.falseNegative_ << std::endl << std::endl;

    std::cout << "testing accuracy: " << gtestAcc << std::endl;
    std::cout << "testing sensitivity: " << testSens << std::endl;
    std::cout << "testing specificity: " << testSpec << std::endl;
    std::cout << "testing precision: " << testPrec << std::endl << std::endl;
    std::cout << "testing true positives: " << gTestQual.truePositive_ << std::endl;
    std::cout << "testing true negatives: " << gTestQual.trueNegative_ << std::endl;
    std::cout << "testing false positives: " << gTestQual.falsePositive_ << std::endl;
    std::cout << "testing false negatives: " << gTestQual.falseNegative_ << std::endl << std::endl;

  }

  std::cout << std::endl;
  std::cout << "===============================================================" << std::endl;
  std::cout << std::endl;

  if (bUseRecursion == false) {
    if (bUseFloat) {
      std::cout << "Needed time: " << gtimings.timeComplete_ << " seconds (Single Precision)" << std::endl;
    } else {
      std::cout << "Needed time: " << gtimings.timeComplete_ << " seconds (Double Precision)" << std::endl;
    }

    std::cout << std::endl << "Timing Details:" << std::endl;
    std::cout << "         mult (complete): " << gtimings.timeMultComplete_ << " seconds" << std::endl;
    std::cout << "         mult (compute) : " << gtimings.timeMultCompute_ << " seconds" << std::endl;
    std::cout << "         mult (comm)    : " << gtimings.timeMultComplete_ - gtimings.timeMultCompute_ << " seconds" << std::endl;
    std::cout << "  mult trans. (complete): " << gtimings.timeMultTransComplete_ << " seconds" << std::endl;
    std::cout << "  mult trans. (compute) : " << gtimings.timeMultTransCompute_ << " seconds" << std::endl;
    std::cout << "  mult trans. (comm)    : " << gtimings.timeMultTransComplete_ - gtimings.timeMultTransCompute_ << " seconds" << std::endl;
    std::cout << std::endl;
    std::cout << "GFlop/s (complete): " << gtimings.GFlop_ / gtimings.timeComplete_ << std::endl;
    std::cout << "GByte/s (complete): " << gtimings.GByte_ / gtimings.timeComplete_ << std::endl;
    std::cout << "GFlop/s (compute): " << gtimings.GFlop_ / (gtimings.timeMultCompute_ + gtimings.timeMultTransCompute_) << std::endl;
    std::cout << "GByte/s (compute): " << gtimings.GByte_ / (gtimings.timeMultCompute_ + gtimings.timeMultTransCompute_) << std::endl << std::endl;
  } else {
    std::cout << "Needed time: " << gtimings.timeComplete_ << " seconds (Double Precision, recursive)" << std::endl << std::endl;
  }

  std::cout << "===============================================================" << std::endl;
  std::cout << std::endl;

  int ompThreadCount = 1;
#ifdef _OPENMP
  #pragma omp parallel
  {
    ompThreadCount = omp_get_num_threads();
  }
#endif
  std::cout << "$" << gdataFile << ";" << gtestFile << ";" << gisRegression << ";" << bUseFloat << ";"
            << ggridtype << ";" << gstart_level << ";" << glambda << ";" << gSLEfinal.maxIterations_ << ";" << gSLEfinal.eps_ << ";"
            << gAdapConfig.numRefinements_ << ";"  << gAdapConfig.threshold_ << ";" << gAdapConfig.noPoints_ << ";"
            << gtrainAcc << ";" << gtestAcc << ";" << gtimings.timeComplete_ << ";" << gtimings.timeMultComplete_
            << ";" << gtimings.timeMultCompute_ << ";" << gtimings.timeMultTransComplete_ << ";" << gtimings.timeMultTransCompute_
            << ";" << gtimings.GFlop_ / gtimings.timeComplete_ << ";"
            << gtimings.GByte_ / gtimings.timeComplete_ << ";" << gtimings.GFlop_ / (gtimings.timeMultCompute_ + gtimings.timeMultTransCompute_) << ";" << gtimings.GByte_ / (gtimings.timeMultCompute_ + gtimings.timeMultTransCompute_)  << ";" << SGPP::parallel::myGlobalMPIComm->getNumRanks()  << ";" << ompThreadCount  << std::endl << std::endl;
}

void adaptClassificationTest(SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataMatrix& testdata, SGPP::base::DataVector& testclasses, bool isRegression,
                             SGPP::base::RegularGridConfiguration& GridConfig, const SGPP::solver::SLESolverConfiguration& SolverConfigRefine,
                             const SGPP::solver::SLESolverConfiguration& SolverConfigFinal, const SGPP::base::AdpativityConfiguration& AdaptConfig,
                             const double lambda, const SGPP::parallel::VectorizationType vecType, const SGPP::parallel::MPIType mpiType) {
  SGPP::datadriven::LearnerBase* myLearner;

  myLearner = new SGPP::parallel::LearnerVectorizedIdentity(vecType, mpiType, isRegression, true);

  // training
  gtimings = myLearner->train(data, classes, GridConfig, SolverConfigRefine,  SolverConfigFinal, AdaptConfig, false, lambda);

  double time_gTrainAcc = 0;
  double time_gTrainQual = 0;
  double time_gTestQual = 0;
  double time_gTestAcc = 0;
  SGPP::base::SGppStopwatch* myStopwatch = new SGPP::base::SGppStopwatch();

  // testing
  myStopwatch->start();
  gtrainAcc = myLearner->getAccuracy(data, classes);
  time_gTrainAcc = myStopwatch->stop();

  myStopwatch->start();
  gtestAcc = myLearner->getAccuracy(testdata, testclasses);
  time_gTestAcc = myStopwatch->stop();

  if (!isRegression) {
    myStopwatch->start();
    gTrainQual = myLearner->getCassificatorQuality(data, classes);
    time_gTrainQual = myStopwatch->stop();

    myStopwatch->start();
    gTestQual = myLearner->getCassificatorQuality(testdata, testclasses);
    time_gTestQual = myStopwatch->stop();
  }

  if (SGPP::parallel::myGlobalMPIComm->getMyRank() == 0) {
    std::cout << "Times for Testing: " << std::endl
              << "Train Acc:  " << time_gTrainAcc << " s" << std::endl
              << "Test  Acc:  " << time_gTestAcc << " s" << std::endl
              << "Train Qual: " << time_gTrainQual << " s" << std::endl
              << "Test  Qual: " << time_gTestQual << " s" << std::endl;
  }

#ifdef GNUPLOT
#endif

  delete myLearner;

  printResults();
}

void adaptClassificationTestRecursive(SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataMatrix& testdata, SGPP::base::DataVector& testclasses, bool isRegression,
                                      SGPP::base::RegularGridConfiguration& GridConfig, const SGPP::solver::SLESolverConfiguration& SolverConfigRefine,
                                      const SGPP::solver::SLESolverConfiguration& SolverConfigFinal, const SGPP::base::AdpativityConfiguration& AdaptConfig,
                                      const double lambda, const SGPP::parallel::VectorizationType vecType) {
  SGPP::datadriven::LearnerBase* myLearner;
  SGPP::datadriven::LearnerRegularizationType C_type;

#ifdef USE_REC_LAPLACE
  C_type = SGPP::datadriven::Laplace;
#else
  C_type = SGPP::datadriven::Identity;
#endif
  myLearner = new SGPP::datadriven::Learner(C_type, isRegression, true);

  // training
  gtimings = myLearner->train(data, classes, GridConfig, SolverConfigRefine,  SolverConfigFinal, AdaptConfig, false, lambda);

  // testing
  gtrainAcc = myLearner->getAccuracy(data, classes);
  gtestAcc = myLearner->getAccuracy(testdata, testclasses);

  if (!isRegression) {
    gTrainQual = myLearner->getCassificatorQuality(data, classes);
    gTestQual = myLearner->getCassificatorQuality(testdata, testclasses);
  }

#ifdef GNUPLOT
#endif

  delete myLearner;

  printResults();
}

void adaptClassificationTestSP(SGPP::base::DataMatrixSP& dataSP, SGPP::base::DataVectorSP& classesSP, SGPP::base::DataMatrixSP& testdataSP, SGPP::base::DataVectorSP& testclassesSP, bool isRegression,
                               SGPP::base::RegularGridConfiguration& GridConfig, const SGPP::solver::SLESolverSPConfiguration& SolverConfigRefine,
                               const SGPP::solver::SLESolverSPConfiguration& SolverConfigFinal, const SGPP::base::AdpativityConfiguration& AdaptConfig,
                               const float lambda, const SGPP::parallel::VectorizationType vecType, const SGPP::parallel::MPIType mpiType) {
  SGPP::datadriven::LearnerBaseSP* myLearner;

  myLearner = new SGPP::parallel::LearnerVectorizedIdentitySP(vecType, mpiType, isRegression, true);

  // training
  gtimings = myLearner->train(dataSP, classesSP, GridConfig, SolverConfigRefine,  SolverConfigFinal, AdaptConfig, false, lambda);

  // testing
  gtrainAcc = myLearner->getAccuracy(dataSP, classesSP);
  gtestAcc = myLearner->getAccuracy(testdataSP, testclassesSP);

  if (!isRegression) {
    gTrainQual = myLearner->getCassificatorQuality(dataSP, classesSP);
    gTestQual = myLearner->getCassificatorQuality(testdataSP, testclassesSP);
  }

#ifdef GNUPLOT
#endif

  delete myLearner;

  printResults();
}

void printHelp() {
  if (SGPP::parallel::myGlobalMPIComm->getMyRank() != 0) {
    return;
  }

  std::cout << std::endl;
  std::cout << "Help for classification/regression benchmark" << std::endl << std::endl;
  std::cout << "Needed parameters:" << std::endl;
  std::cout << "	Traindata-file" << std::endl;
  std::cout << "	Testdata-file" << std::endl;
  std::cout << "	regression (0/1)" << std::endl;
  std::cout << "	precision (SP,DP)" << std::endl;
  std::cout << "	gridtype (linear,linearboundary,modlinear)" << std::endl;
  std::cout << "	Startlevel" << std::endl;
  std::cout << "	lambda" << std::endl;
  std::cout << "	CG max. iterations" << std::endl;
  std::cout << "	CG epsilon" << std::endl;
  std::cout << "	#refinements" << std::endl;
  std::cout << "	Refinement threshold" << std::endl;
  std::cout << "	#points refined" << std::endl;
  std::cout << "	CG max. iterations, first refinement steps" << std::endl;
  std::cout << "	CG epsilon, first refinement steps" << std::endl;
  std::cout << "	Vectorization: X86SIMD, OCL, HYBRID_X86SIMD_OCL, ArBB; " << std::endl;
  std::cout << "			for classical sparse grid algorithms choose: REC" << std::endl << std::endl << std::endl;
  std::cout << "	MPI Communication Method: NONE, Allreduce, Alltoallv, Async, Onesided, TrueAsync, Bigdata; " << std::endl;
  std::cout << "Example call:" << std::endl;
  std::cout << "	app.exe     test.data train.data 0 SP linearboundary 3 0.000001 250 0.0001 6 0.0 100 20 0.1 X86SIMD" << std::endl << std::endl << std::endl;
}

void printHeader() {
  if (SGPP::parallel::myGlobalMPIComm->getMyRank() != 0) {
    return;
  }

  std::cout << std::endl;
  std::cout << "===============================================================" << std::endl;
  std::cout << "Classification Test Application" << std::endl;
  std::cout << "===============================================================" << std::endl;
  std::cout << std::endl;
}

/**
 * Testapplication for the Intel VTune Profiling Tool
 * and a measurement app for Sparse Grid Algorithms building blocks
 */
int main(int argc, char* argv[]) {
  int mpi_myid;
  int mpi_size;

  int threadLevelProvided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &threadLevelProvided);

  if (threadLevelProvided != MPI_THREAD_FUNNELED) {
    std::cout << "MPI Library does not support Multithreaded Processes" << MPI_Finalize();
    return -1;
  }

  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myid);
  SGPP::parallel::myGlobalMPIComm = new SGPP::parallel::MPICommunicator(mpi_myid, mpi_size);

  std::streambuf* stdoutBuf = std::cout.rdbuf();
  std::ofstream dummy_out("/dev/null");

  if (mpi_myid != 0) { // disable output for all processes but proc 0
    std::cout.rdbuf(dummy_out.rdbuf());
  }

  //std::cout << "Startup of Process " << mpi_myid << std::endl;

  std::string dataFile;
  std::string testFile;
  std::string gridtype;
  std::string precision;
  std::string vectorization;
  std::string mpiConfValue;

  double lambda = 0.0;
  double cg_eps = 0.0;
  double refine_thresh = 0.0;
  double cg_eps_learning = 0.0;

  size_t cg_max = 0;
  size_t refine_count = 0;
  size_t refine_points = 0;
  int start_level = 0;
  size_t cg_max_learning = 0;

  SGPP::base::RegularGridConfiguration gridConfig;
  SGPP::solver::SLESolverConfiguration SLESolverConfigRefine;
  SGPP::solver::SLESolverConfiguration SLESolverConfigFinal;
  SGPP::solver::SLESolverSPConfiguration SLESolverSPConfigRefine;
  SGPP::solver::SLESolverSPConfiguration SLESolverSPConfigFinal;
  SGPP::base::AdpativityConfiguration adaptConfig;
  SGPP::parallel::VectorizationType vecType;
  SGPP::parallel::MPIType mpiType;

  bool regression;

  printHeader();

  if (argc != 17) {
    printHelp();
  } else {
    dataFile.assign(argv[1]);
    testFile.assign(argv[2]);
    regression = false;

    if (atoi(argv[3]) == 1) {
      regression = true;
    }

    precision.assign(argv[4]);
    gridtype.assign(argv[5]);
    start_level = atoi(argv[6]);
    lambda = atof(argv[7]);
    cg_max = atoi(argv[8]);
    cg_eps = atof(argv[9]);
    refine_count = atoi(argv[10]);
    refine_thresh = atof(argv[11]);
    refine_points = atoi(argv[12]);
    cg_max_learning = atoi(argv[13]);
    cg_eps_learning = atof(argv[14]);
    vectorization.assign(argv[15]);
    mpiConfValue.assign(argv[16]);

    // Set Vectorization
    // Fallback
    if (vectorization == "X86SIMD") {
      vecType = SGPP::parallel::X86SIMD;
    } else if  (vectorization == "OCL") {
      vecType = SGPP::parallel::OpenCL;
    } else if  (vectorization == "HYBRID_X86SIMD_OCL") {
      vecType = SGPP::parallel::Hybrid_X86SIMD_OpenCL;
    } else if  (vectorization == "ArBB") {
      vecType = SGPP::parallel::ArBB;
    } else if  (vectorization == "MIC") {
      vecType = SGPP::parallel::MIC;
    } else if  (vectorization == "HYBRID_X86SIMD_MIC") {
      vecType = SGPP::parallel::Hybrid_X86SIMD_MIC;
    } else {
      vecType = SGPP::parallel::X86SIMD;
    }

    // set MPI Type
    if (mpiConfValue == "Alltoallv") {
      mpiType = SGPP::parallel::MPIAlltoallv;
    } else if (mpiConfValue == "Allreduce") {
      mpiType = SGPP::parallel::MPIAllreduce;
    } else if (mpiConfValue == "Async") {
      mpiType = SGPP::parallel::MPIAsync;
    } else if (mpiConfValue == "TrueAsync") {
      mpiType = SGPP::parallel::MPITrueAsync;
    } else if (mpiConfValue == "Onesided") {
      mpiType = SGPP::parallel::MPIOnesided;
    } else if (mpiConfValue == "Bigdata") {
      mpiType = SGPP::parallel::MPIBigdata;
    } else {
      mpiType = SGPP::parallel::MPINone;
    }

    // Set Adaptivity
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = refine_points;
    adaptConfig.numRefinements_ = refine_count;
    adaptConfig.percent_ = 100.0;
    adaptConfig.threshold_ = refine_thresh;

    // Set solver during refinement
    SLESolverConfigRefine.eps_ = cg_eps_learning;
    SLESolverConfigRefine.maxIterations_ = cg_max_learning;
    SLESolverConfigRefine.threshold_ = -1.0;
    SLESolverConfigRefine.type_ = SGPP::solver::CG;

    SLESolverSPConfigRefine.eps_ = static_cast<float>(cg_eps_learning);
    SLESolverSPConfigRefine.maxIterations_ = cg_max_learning;
    SLESolverSPConfigRefine.threshold_ = -1.0f;
    SLESolverSPConfigRefine.type_ = SGPP::solver::CG;

    // Set solver for final step
    SLESolverConfigFinal.eps_ = cg_eps;
    SLESolverConfigFinal.maxIterations_ = cg_max;
    SLESolverConfigFinal.threshold_ = -1.0;
    SLESolverConfigFinal.type_ = SGPP::solver::CG;

    SLESolverSPConfigFinal.eps_ = static_cast<float>(cg_eps);
    SLESolverSPConfigFinal.maxIterations_ = cg_max;
    SLESolverSPConfigFinal.threshold_ = -1.0f;
    SLESolverSPConfigFinal.type_ = SGPP::solver::CG;

    std::string tfileTrain = dataFile;
    std::string tfileTest = testFile;

    SGPP::datadriven::Dataset dataset;
    SGPP::datadriven::Dataset testdataset = SGPP::datadriven::ARFFTools::readARFF(tfileTest);

    size_t nDim;
    size_t nInstancesNo;
    size_t nInstancesTestNo;

    SGPP::datadriven::DatasetGenerator* g = NULL;

    if (mpiType == SGPP::parallel::MPIBigdata) {
      if (dataFile.find("friedman1") != std::string::npos) {
        g = new SGPP::datadriven::Friedman1Generator();
      } else if (dataFile.find("friedman2") != std::string::npos) {
        g = new SGPP::datadriven::Friedman2Generator();
      } else if (dataFile.find("friedman3") != std::string::npos) {
        g = new SGPP::datadriven::Friedman3Generator();
      } else {
        std::cout << "cannot generate dataset for " << dataFile << std::endl;
        throw SGPP::base::operation_exception("cannot generate dataset");
      }

      nDim = g->getDims();
      nInstancesNo = 100000; // number of instances per node

      const char* dataset_generation_count = getenv("SGPP_DATASET_GENERATION_COUNT");

      if (dataset_generation_count != NULL) {
        nInstancesNo = (size_t)(strtoul (dataset_generation_count, NULL, 0));
      }

      std::cout << "Generating " << nInstancesNo << " datasets per node (for " << mpi_size << " nodes)." << std::endl;

    } else {
      SGPP::datadriven::ARFFTools::readARFFSize(tfileTrain, nInstancesNo, nDim);
    }

    nInstancesTestNo = testdataset.getNumberInstances();

    // Define DP data
    SGPP::base::DataMatrix data;
    SGPP::base::DataVector classes;
    SGPP::base::DataMatrix testdata = *testdataset.getTrainingData();
    SGPP::base::DataVector testclasses = *testdataset.getClasses();

    // Define SP data
    SGPP::base::DataMatrixSP dataSP(nInstancesNo, nDim);
    SGPP::base::DataVectorSP classesSP(nInstancesNo);
    SGPP::base::DataMatrixSP testdataSP(nInstancesTestNo, nDim);
    SGPP::base::DataVectorSP testclassesSP(nInstancesTestNo);


    if (mpiType == SGPP::parallel::MPIBigdata) {
      g->createData(mpi_myid, nInstancesNo, data, classes);
      delete g;
    } else {
      // Read data from file
      dataset = SGPP::datadriven::ARFFTools::readARFF(tfileTrain);
      data = *dataset.getTrainingData();
      classes = *dataset.getClasses();
    }

    SGPP::base::PrecisionConverter::convertDataMatrixToDataMatrixSP(data, dataSP);
    SGPP::base::PrecisionConverter::convertDataVectorToDataVectorSP(classes, classesSP);
    SGPP::base::PrecisionConverter::convertDataMatrixToDataMatrixSP(testdata, testdataSP);
    SGPP::base::PrecisionConverter::convertDataVectorToDataVectorSP(testclasses, testclassesSP);

    // Set Grid-Information
    gridConfig.dim_ = nDim;
    ggridtype = gridtype;

    if (gridtype == "linearboundary") {
      gridConfig.type_ = SGPP::base::LinearBoundary;
    } else if (gridtype == "modlinear") {
      gridConfig.type_ = SGPP::base::ModLinear;
    } else if (gridtype == "linear") {
      gridConfig.type_ = SGPP::base::Linear;
    } else {
      std::cout << std::endl << "An unsupported grid type was chosen! Exiting...." << std::endl << std::endl;
      return -1;
    }

    gridConfig.level_ = start_level;
    glambda = lambda;
    gisRegression = regression;
    gdataFile = dataFile;
    gtestFile = testFile;
    gSLEfinal = SLESolverConfigFinal;
    gAdapConfig = adaptConfig;
    gstart_level = start_level;


    if (mpi_myid == 0) {
      std::cout << std::endl << "Dims: " << nDim << "; Traininstances: " << nInstancesNo << "; Testinstances: " << nInstancesTestNo << std::endl << std::endl;
    }

    if (vectorization == "REC") {
      bUseFloat = false;
      bUseRecursion = true;

      printSettings(dataFile, testFile, regression, gridConfig, SLESolverConfigRefine,
                    SLESolverConfigFinal, adaptConfig, lambda, vecType);

      adaptClassificationTestRecursive(data, classes, testdata, testclasses, regression, gridConfig, SLESolverConfigRefine,
                                       SLESolverConfigFinal, adaptConfig, lambda, vecType);
    } else {
      bUseRecursion = false;

      if (precision == "SP") {
        bUseFloat = true;

        printSettings(dataFile, testFile, regression, gridConfig, SLESolverConfigRefine,
                      SLESolverConfigFinal, adaptConfig, lambda, vecType);

        adaptClassificationTestSP(dataSP, classesSP, testdataSP, testclassesSP, regression, gridConfig, SLESolverSPConfigRefine,
                                  SLESolverSPConfigFinal, adaptConfig, (float)lambda, vecType, mpiType);
      } else if (precision == "DP") {
        bUseFloat = false;

        printSettings(dataFile, testFile, regression, gridConfig, SLESolverConfigRefine,
                      SLESolverConfigFinal, adaptConfig, lambda, vecType);

        adaptClassificationTest(data, classes, testdata, testclasses, regression, gridConfig, SLESolverConfigRefine,
                                SLESolverConfigFinal, adaptConfig, lambda, vecType, mpiType);
      } else {
        std::cout << "Unsupported precision type has been chosen! Existing...." << std::endl << std::endl;
        return -1;
      }
    }
  }

  if (mpi_myid != 0) { // restore stdout buffer
    std::cout.rdbuf(stdoutBuf);
  }

  MPI_Finalize();

  return 0;
}
