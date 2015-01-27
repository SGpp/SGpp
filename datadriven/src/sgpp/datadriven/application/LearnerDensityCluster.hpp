// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef LEARNERDENSITYCLUSTER_HPP
#define LEARNERDENSITYCLUSTER_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/application/LearnerBase.hpp>
#include <sgpp/datadriven/tools/TypesDatadriven.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <vector>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {
	/**
	 * enum to select different threshold function
	 */
	enum ThresholdType {
      Constant,
      Relative,
      Difference,
      Minima
    };

    /**
     * structure that can be used by applications to cluster density-based clustering information
     */
    struct DensityBasedClusteringConfiguration {
	  /// Number of neighbors
	  int numberOfNeighbors;
	  /// Epsilon
	  double eps;
	  /// Method that is used to remove edges from the graph
	  SGPP::datadriven::ThresholdType thresholdType;
	  /// Threshold value
	  double threshold;
    };


    class LearnerDensityCluster: public SGPP::datadriven::LearnerBase {
      protected:
    	/// contains the clustering result
    	SGPP::base::DataVector* clusterAssignments_ = NULL;
    	/// contains the evaluated grid values for the data points
    	SGPP::base::DataVector* gridVals_ = NULL;
    	/// contains the calculated nearest neighbors
    	int ** neighbors_ = NULL;

    	int neighborsRows_ = 0;
    	int neighborsCols_ = 0;

    	/// epsilon for the threshold function
    	double eps = 5;
    	/// threshold for the local minima
    	double threshold = 0;
    	/// how many neighbors should taken into consideration
    	int numberOfNeighbors = 5;
    	/// specify the threshold function
    	SGPP::datadriven::ThresholdType thresholdType = SGPP::datadriven::Constant;

    	/// it is used to save whether a point is a minimum point or not
    	std::vector<bool> * minimumPoint_ = NULL;

    	/**
    	 * evalute the grid for the given points
    	 *
    	 * @param testDataset the given data
    	 * @param GridConfig the grid configuration
    	 * @param SolverConfig the solver configuration
    	 * @param lambda the regularization parameter
    	 */
    	void calculateGridValues(SGPP::base::DataMatrix& testDataset, const SGPP::base::RegularGridConfiguration& GridConfig,
    			const SGPP::solver::SLESolverConfiguration& SolverConfig, const double lambda);

    	/**
    	 * save a double array
    	 *
    	 * @param tFilename where to save the array
    	 * @param len the length of the array
    	 * @param data the array
    	 */
    	void saveArray(const char * tFilename, int len, double data[]);

    	/**
    	 * load a double array
    	 *
    	 * @param tFilename from where to load
    	 * @param len the length of the array
    	 * @return the array from the file
    	 */
    	double * loadArray(const char * tFilename,int * len);

    	/**
    	 * save a two-dimensional double array
    	 *
    	 * @param tFilename where to save the array
    	 * @param row how many rows
    	 * @param col how many columns
    	 * @param data the two-dimensional array
    	 */
    	void saveArray(const char * tFilename, int row, int col, int ** data);

    	/**
    	 * load a two-dimensional double array
    	 *
    	 * @param tFilename from where to load
    	 * @param row how many rows
    	 * @param col how many columns
    	 * @return the two-dimensional array
    	 */
    	int ** loadArray(const char * tFilename, int * row, int * col);

    	/**
    	 * constant threshold function (f_t1)
    	 *
    	 * @param testDataset the data
    	 * @param i the index of the first point
    	 * @param j the index of the second point
    	 * @return boolean whether the points should be connected or not
    	 */
    	bool constantThreshold(SGPP::base::DataMatrix& testDataset, int i, int j);

    	/**
    	 * relative threshold function (f_t2)
    	 *
    	 * @param testDataset the data
    	 * @param i the index of the first point
    	 * @param j the index of the second point
    	 * @return boolean whether the points should be connected or no
    	 */
    	bool relativeThreshold(SGPP::base::DataMatrix& testDataset, int i, int j);

    	/**
		 * difference threshold function (f_t2)
		 *
		 * @param testDataset the data
		 * @param i the index of the first point
		 * @param j the index of the second point
		 * @return boolean whether the points should be connected or no
		 */
    	bool differenceThreshold(SGPP::base::DataMatrix& testDataset, int i, int j);

    	/**
    	 * delete the neighbors in neighbors_
    	 */
    	void deleteNeighbors();

    	/**
    	 * compare two pairs
    	 *
    	 * @param firstElem
    	 * @param secondElem
    	 * @return
    	 */
    	static bool pairCompare(const std::pair<int, double>& firstElem, const std::pair<int, double>& secondElem);
      public:

    	/**
    	 * Default constructor
    	 */
    	LearnerDensityCluster();

    	/**
    	 * Constructor
    	 *
    	 * @param isVerbose
    	 */
    	LearnerDensityCluster(bool isVerbose);

    	/**
    	 * Destructor
    	 */
    	~LearnerDensityCluster();

<<<<<<< .mine
    	/**
    	 * The cluster function. Before the clustering the cluster configuration must be set with the function ::setClusterConfiguration
    	 *
    	 * @param testDataset the given data
    	 * @param classes not used
    	 * @param GridConfig the grid configuration
    	 * @param SolverConfig the solver configuration
    	 * @param lambda the regularization parameter
    	 * @return TODO: empty
    	 */
    	SGPP::datadriven::LearnerTiming train(SGPP::base::DataMatrix& testDataset, SGPP::base::DataVector& classes,
    						  const SGPP::base::RegularGridConfiguration& GridConfig, const SGPP::solver::SLESolverConfiguration& SolverConfig,
    						  const double lambda);
=======
        /**
         * The cluster function. Before the clustering the cluster configuration must be set with the function LearnerDensityCluster::setClusterConfiguration
         *
         * @param testDataset the given data
         * @param classes not used
         * @param GridConfig the grid configuration
         * @param SolverConfig the solver configuration
         * @param lambda the regularization parameter
         * @return TODO: empty
         */
        SGPP::datadriven::LearnerTiming train(SGPP::base::DataMatrix& testDataset, SGPP::base::DataVector& classes,
                                              const SGPP::base::RegularGridConfiguration& GridConfig, const SGPP::solver::SLESolverConfiguration& SolverConfig,
                                              const double lambda);
>>>>>>> .r4048

    	/**
    	 * TODO: Not implemented
    	 * @param testDataset
    	 * @param lambda
    	 * @return
    	 */
    	SGPP::datadriven::DMSystemMatrixBase* createDMSystem(SGPP::base::DataMatrix& testDataset, double lambda);

    	/**
    	 * Gets the last clustering result
    	 *
    	 * @return clustering result
    	 */
    	SGPP::base::DataVector* getClusterAssignments();

    	/**
    	 * Precalculates the grid value for each point and stores the information a file.
    	 *
    	 * @param filename where to store
    	 * @param testDataset the data
    	 * @param GridConfig the grid configuration
    	 * @param SolverConfig the solver configuration
    	 * @param lamda the regularization parameter
    	 */
    	void precalculateGridValues(const char * filename, SGPP::base::DataMatrix& testDataset, const SGPP::base::RegularGridConfiguration& GridConfig, const SGPP::solver::SLESolverConfiguration& SolverConfig,
				  const double lamda);

    	/**
    	 * Load the precalculated values from a file
    	 * @param filename the file
    	 */
    	void loadPrecalculatedValues(const char * filename);

<<<<<<< .mine
    	/**
    	 * Cluster function. The grid values must be calculated with ::precalculateGridValues and the cluster configuration must be set with the function ::setClusterConfiguration
    	 * @param trainDataset the data
    	 * @param GridConfig the grid configuration
    	 */
    	void cluster(SGPP::base::DataMatrix& trainDataset, const SGPP::base::RegularGridConfiguration& GridConfig);
=======
        /**
         * Cluster function. The grid values must be calculated with LearnerDensityCluster::precalculateGridValues and the cluster configuration must be set with the function LearnerDensityCluster::setClusterConfiguration
         * @param trainDataset the data
         * @param GridConfig the grid configuration
         */
        void cluster(SGPP::base::DataMatrix& trainDataset, const SGPP::base::RegularGridConfiguration& GridConfig);
>>>>>>> .r4048

    	/**
    	 * Set the cluster configuration
    	 *
    	 * @param ClusterConfig the cluster configuration
    	 */
    	void setClusterConfiguration(const SGPP::datadriven::DensityBasedClusteringConfiguration& ClusterConfig);

    	/**
    	 * Do a simple classification with the clustering result.The n biggest clusters form the classes.
    	 *
    	 * @param trainDataset the data set
    	 * @param components the clustering result
    	 * @param n
    	 * @param newComponents
    	 */
        void postprocessing(SGPP::base::DataMatrix& trainDataset,SGPP::base::DataVector& components, int n, SGPP::base::DataVector& newComponents);

    	/**
    	 * Prcalculates the n nearest neighbors and store the information in a file
    	 * @param filename the file
    	 * @param testDataset the data
    	 * @param n n nearest neighbors
    	 */
    	void precalculateNeighbors(const char * filename, SGPP::base::DataMatrix& testDataset, int n);

    	/**
    	 * Loads the precalculated neighbors
    	 * @param filename file
    	 */
    	void loadPrecalculatedNeighbors(const char * filename);
    };

  }

}

#endif /* LEARNERDENSITYCLUSTER_HPP */
