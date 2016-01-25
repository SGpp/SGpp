// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/solver/SLESolver.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>

#include <deque>
#include <map>
#include "BatchConfiguration.hpp"

using namespace SGPP::base;
using namespace std;

namespace SGPP {
  namespace datadriven {
    /**
     *  The Batchlearner learns the data provided as input in batches. The batches can by weighted with different functions. Labels are predicted by using density functions. Adaptivity can be enabled.
     */
    class BatchLearner {
      protected:
        SGPP::base::BatchConfiguration batchConf;//!< configuration for the BatchLearner
        SGPP::base::RegularGridConfiguration gridConf;//!< configuration for the grids
        SGPP::solver::SLESolverConfiguration solverConf;//!< configuration for the solver
        SGPP::base::AdpativityConfiguration adaptConf;//!< configuration for the adaptivity
        std::fstream reader;//!< stream to read in the arff file
        bool reachedData = false;//!< flag if "\@DATA" has been reached in the arff
        int batchnum = 0;//!< number of the current batch
        size_t dimensions = 0;//!< count of dimensions in the data
        std::map<int, SGPP::base::DataMatrix*> dataInBatch;//!< mapping of data in each batch to label
        std::map<int, SGPP::base::DataVector*> alphaVectors;//!< mapping of alpha vectors to label
        std::map<int, SGPP::base::LinearGrid*> grids;//!< mapping of grids to label
        std::map<int, float_t> normFactors;//!< mapping of factors for the normalization to label
        std::map<int, size_t> occurences;//!< mapping of count of items to label
        std::map<int, std::deque<SGPP::base::DataVector> > alphaStorage;//!< mapping used to store the previous alphas

        int dataLine = 0;//!< stores the current line number from \@DATA
        unsigned int bs = 0; //!< count of lines read for the batch
        std::string batch; //!< because testdata is added to the following batch: save batch global
        int t_total = 0;//!< total items tested
        int t_correct = 0;//!< items predicted correct
        SGPP::solver::SLESolver* myCG;//!< solver

        /**
         * function to apply one of the weighting methods to the alpha provided, returns the weighted alpha
         * @see batchConfig.wMode for the different weighting methods
         * @param alpha recent alpha that has to be taken into account for the calculation
         * @param grid the class the alphas referes to
         */
        SGPP::base::DataVector applyWeight(SGPP::base::DataVector alpha, int grid);

        /**
         * function that processes the data provided as string, learns the new data and refines the grids (if wanted by user)
         * @param workData string containing the items separated by '\n', data is separated by ',', the last entry per item is the class as int
         */
        void processBatch(std::string workData);

        /**
         * function that parses a string to data and class
         * @param input string of one item, values separated by ',', last entry is the class as int
         * @param dataFound the data found in the string
         * @param classFound the class found in the string
         */
        void stringToDataVector(std::string input, SGPP::base::DataVector& dataFound, int& classFound);

        /**
         * function that parses a string containing many items
         * if (mapData) and maps the data according to class into dataInBatch (will be cleared before usage, dataFound and classesFound will not be touched)
         * if (!mapData) to a DataMatrix containg the data and a DataVector containing the classes (dataInBatch will not be touched)
         * @param input many items seperated by '\n'
         * @param dataFound DataMatrix that will contain the found data afterwards (will be cleared before usage)
         * @param classesFound DataVector that will contain the found classes afterwards (will be cleared before usage)
         * @param mapData controls whether data is mapped to dataInBatch or saved to dataFound and classesFound
        */
        void stringToDataMatrix(std::string& input, SGPP::base::DataMatrix& dataFound, SGPP::base::DataVector& classesFound, bool mapData);

        bool isFinished = false; //!< indicates whether the stream has been read to the end
        float_t acc_global = -1.0; //!< accuracy over all predictions done so far (including the ones if batchConfig.testsize > 0 )
        float_t acc_current = -1.0; //!< accuracy of the last call of predict(..)
      public:
        /**
         * constructor taking all relevant parameters
         * @param batchConfig struct containig all parameters specific for the BatchLearner
         * @param gridConfig configuration for the grids
         * @param solverConfig configuration for the solver
         * @param adaptConfig configuration for the adaptivity of the solver
         */
        BatchLearner(
          SGPP::base::BatchConfiguration batchConfig,
          SGPP::base::RegularGridConfiguration gridConfig,
          SGPP::solver::SLESolverConfiguration solverConfig,
          SGPP::base::AdpativityConfiguration adaptConfig);

        /**
         * close the stream
         */
        void closeFile() {
          reader.close();
        }

        /**
         * Get whether stream has been read to the end
         */
        bool getIsFinished() {
          return isFinished;
        }

        /**
                           * Get the accuracy of the last batch predicted
                           */
        float_t getAccCurrent() {
          return acc_current;
        }

        /**
                           * Get the accuracy over all predictions
                           */
        float_t getAccGlobal() {
          return acc_global;
        }

        /**
         * predict labels of the data provided, return vector of predicted labels
         * @param entries DataMatrix containing the data to test
         * @param updateNorm should the normalization factors be updated before predicting?
               */
        SGPP::base::DataVector predict(SGPP::base::DataMatrix& entries, bool updateNorm);

        /**
               * read all lines needed for one batch (and maybe test data), call processBatch(..)
         * if test is wanted by user: call predict(..) on test data, calculate, save and output accuracies
               */
        void trainBatch();

    };
  }
}
