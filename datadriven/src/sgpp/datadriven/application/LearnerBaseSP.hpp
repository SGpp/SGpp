// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef LEARNERBASESP_HPP
#define LEARNERBASESP_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVectorSP.hpp>
#include <sgpp/base/datatypes/DataMatrixSP.hpp>
#include <sgpp/base/tools/PrecisionConverter.hpp>

#include <sgpp/solver/SLESolverSP.hpp>

#include <sgpp/datadriven/algorithm/DMSystemMatrixBaseSP.hpp>
#include <sgpp/datadriven/tools/TypesDatadriven.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace datadriven {

    /**
     * Abstract class that implements a regression/classification learner
     * based on spatial adaptive Sparse Grids.
     *
     * Furthermore this class is intended to provide a framework
     * for adavanded regression and classification methods
     * by allowing to override basic routines like train
     * and test.
     *
     * This versions supports single precision datatypes.
     */
    class LearnerBaseSP {
      protected:
        /// the grid's coefficients
        SGPP::base::DataVectorSP* alpha_;
        /// sparse grid object
        SGPP::base::Grid* grid_;
        /// is verbose output enabled
        bool isVerbose_;
        /// is regression selected
        bool isRegression_;
        /// is the grid trained
        bool isTrained_;
        /// execution time
        double execTime_;
        /// number of executed Floating Point operations
        double GFlop_;
        /// number of transferred Gbytes
        double GByte_;

        /**
         * Hook-Method for pre-processing before
         * starting learning.
         *
         * can be overwritten by derived classes
         */
        virtual void preProcessing();

        /**
         * Hook-Method for post-processing after each
         * refinement learning.
         *
         * can be overwritten by derived classes
         * @param trainDataset matrix with training data set
         * @param solver solver
         * @param numNeededIterations number of required iterations
         *
         */
        virtual void postProcessing(const SGPP::base::DataMatrixSP& trainDataset, const SGPP::solver::SLESolverType& solver,
                                    const size_t numNeededIterations);

        /**
         * Initialize the grid and its coefficients
         *
         * @param GridConfig structure which describes the regular start grid
         */
        virtual void InitializeGrid(const SGPP::base::RegularGridConfiguration& GridConfig);

        /**
         * abstract method that constructs the corresponding system of linear equations
         * Derived classes MUST overwrite this functions!
         *
         * @param trainDataset training dataset
         * @param lambda lambda regularization parameter
         */
        virtual SGPP::datadriven::DMSystemMatrixBaseSP* createDMSystem(SGPP::base::DataMatrixSP& trainDataset, float lambda) = 0;

      public:
        /**
         * Constructor
         *
         * @param isRegression flag for regression
         * @param isVerbose flag for verbose output
         */
        LearnerBaseSP(const bool isRegression, const bool isVerbose = true);

        /**
         * Constructor
         *
         * @param tGridFilename path to file that contains a serialized grid
         * @param tAlphaFilename path to file that contains the grid's coefficients
         * @param isRegression set to true if a regression task should be executed
         * @param isVerbose set to true in order to allow console output
         */
        LearnerBaseSP(std::string tGridFilename, std::string tAlphaFilename, const bool isRegression, const bool isVerbose = true);

        /**
         * Copy-Constructor
         *
         * @param copyMe LearnerBase that should be duplicated
         */
        LearnerBaseSP(const LearnerBaseSP& copyMe);

        /**
         * Destructor
         */
        virtual ~LearnerBaseSP();

        /**
         * Learning a dataset with spatially adaptive sparse grids
         *
         * @param testDataset the training dataset
         * @param classes classes corresponding to the training dataset
         * @param GridConfig configuration of the regular start grid
         * @param SolverConfigRefine configuration of the SLE solver during the adaptive refinements of the grid
         * @param SolverConfigFinal configuration of the final SLE solving step on the refined grid
         * @param AdaptConfig configuration of the adaptivity strategy
         * @param testAccDuringAdapt set to true if the training accuracy should be determined in evert refinement step
         * @param lambda regularization parameter lambda
         */
        virtual LearnerTiming train(SGPP::base::DataMatrixSP& testDataset, SGPP::base::DataVectorSP& classes,
                                    const SGPP::base::RegularGridConfiguration& GridConfig, const SGPP::solver::SLESolverSPConfiguration& SolverConfigRefine,
                                    const SGPP::solver::SLESolverSPConfiguration& SolverConfigFinal, const SGPP::base::AdpativityConfiguration& AdaptConfig,
                                    bool testAccDuringAdapt, const float lambda);

        /**
         * Learning a dataset with regular sparse grids
         *
         * @param testDataset the training dataset
         * @param classes classes corresponding to the training dataset
         * @param GridConfig configuration of the regular grid
         * @param SolverConfig configuration of the SLE solver
         * @param lamda regularization parameter lambda
         */
        LearnerTiming train(SGPP::base::DataMatrixSP& testDataset,  SGPP::base::DataVectorSP& classes,
                            const SGPP::base::RegularGridConfiguration& GridConfig, const SGPP::solver::SLESolverSPConfiguration& SolverConfig,
                            const float lamda);

        /**
         * executes a Regression test for a given dataset and returns the result
         *
         * @param testDataset dataset that is evaluated with the current learner
         * @return regression values of testDataset
         */
        virtual SGPP::base::DataVectorSP predict(SGPP::base::DataMatrixSP& testDataset);

        /**
         * compute the accuracy for given testDataset. test is automatically called
         * in order to determine the regression values of the current learner
         *
         * In case if classification (isRegression == false) this routine returns the learner's accuracy
         * In case of regressions (isRegression == true) this routine returns the learner's MSE
         *
         * @param testDataset dataset to be tested
         * @param classesReference reference labels of the test dataset
         * @param threshold threshold used for classification, ignored when performing regressions
         *
         * @return accuracy, percent or MSE, depending on the execution mode
         */
        virtual double getAccuracy(SGPP::base::DataMatrixSP& testDataset,
                                   const SGPP::base::DataVectorSP& classesReference, const float threshold = 0.0);

        /**
         * compute the accuracy for given testDataset.
         *
         * In case if classification (isRegression == false) this routine returns the learner's accuracy
         * In case of regressions (isRegression == true) this routine returns the learner's MSE
         *
         * @param classesComputed regression results of the test dataset
         * @param classesReference reference labels of the test dataset
         * @param threshold threshold used for classification, ignored when performing regressions
         *
         * @return accuracy, percent or MSE, depending on the execution mode
         */
        virtual double getAccuracy(const SGPP::base::DataVectorSP& classesComputed,
                                   const SGPP::base::DataVectorSP& classesReference, const float threshold = 0.0);

        /**
         * compute the quality for given testDataset, classification ONLY!
         * test is automatically called
         * in order to determine the regression values of the current learner
         *
         * @param testDataset dataset to be tested
         * @param classesReference reference labels of the test dataset
         * @param threshold threshold used for classification, ignored when performing regressions
         *
         * @return quality structure containing tp, tn, fp, fn counts
         */
        virtual ClassificatorQuality getCassificatorQuality(SGPP::base::DataMatrixSP& testDataset,
            const SGPP::base::DataVectorSP& classesReference, const float threshold = 0.0);

        /**
         * compute the quality for given testDataset, classification ONLY!
         *
         * @param classesComputed regression results of the test dataset
         * @param classesReference reference labels of the test dataset
         * @param threshold threshold used for classification, ignored when performing regressions
         *
         * @return quality structure containing tp, tn, fp, fn counts
         */
        virtual ClassificatorQuality getCassificatorQuality(const SGPP::base::DataVectorSP& classesComputed,
            const SGPP::base::DataVectorSP& classesReference, const float threshold = 0.0);

        /**
         * store the grid and its current coefficients into files for
         * further usage.
         *
         * @param tGridFilename filename of grid file
         * @param tAlphaFilename filename of coefficient file
         */
        void store(std::string tGridFilename, std::string tAlphaFilename);

        /**
         * simple dump of grid points into file, e.g. used to
         * plot with gnuplot
         *
         * only executed if grid is trained
         *
         * @param tFilename filename to store the dump to
         */
        void dumpGrid(std::string tFilename);

        /**
         * simple dump of sparse grid function into file, e.g. used to
         * plot with gnuplot.
         *
         * only executed if grid is trained and number of dimensions <= 2.
         *
         * @param tFilename filename to store the dump to
         * @param resolution resolution of function plot
         */
        void dumpFunction(std::string tFilename, size_t resolution);

        /**
         * determines the current mode
         *
         * @return returns whether the current mode is regression or not
         */
        bool getIsRegression() const;

        /**
         * determines the current verbose mode of learner
         *
         * @return returns whether the current learner has verbose output
         */
        bool getIsVerbose() const;

        /**
         * sets the current verbose mode of learner
         *
         * @param isVerbose the current learner's verbose output
         */
        void setIsVerbose(const bool isVerbose);
    };

  }

}

#endif /* LEARNERBASESP_HPP */