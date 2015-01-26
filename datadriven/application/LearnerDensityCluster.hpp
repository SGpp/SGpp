/******************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Florian Zipperle (florian.zipperle@tum.de)

#ifndef LEARNERDENSITYCLUSTER_HPP
#define LEARNERDENSITYCLUSTER_HPP

#include "base/datatypes/DataVector.hpp"
#include "datadriven/application/LearnerBase.hpp"
#include "datadriven/tools/TypesDatadriven.hpp"
#include "datadriven/algorithm/DensitySystemMatrix.hpp"
#include <vector>

namespace sg {
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
	  sg::datadriven::ThresholdType thresholdType;
	  /// Threshold value
	  double threshold;
    };


    class LearnerDensityCluster: public sg::datadriven::LearnerBase {
      protected:
    	/// contains the clustering result
    	sg::base::DataVector* clusterAssignments_ = NULL;
    	/// contains the evaluated grid values for the data points
    	sg::base::DataVector* gridVals_ = NULL;
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
    	sg::datadriven::ThresholdType thresholdType = sg::datadriven::Constant;

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
    	void calculateGridValues(sg::base::DataMatrix& testDataset, const sg::base::RegularGridConfiguration& GridConfig,
    			const sg::solver::SLESolverConfiguration& SolverConfig, const double lambda);

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
    	bool constantThreshold(sg::base::DataMatrix& testDataset, int i, int j);

    	/**
    	 * relative threshold function (f_t2)
    	 *
    	 * @param testDataset the data
    	 * @param i the index of the first point
    	 * @param j the index of the second point
    	 * @return boolean whether the points should be connected or no
    	 */
    	bool relativeThreshold(sg::base::DataMatrix& testDataset, int i, int j);

    	/**
		 * difference threshold function (f_t2)
		 *
		 * @param testDataset the data
		 * @param i the index of the first point
		 * @param j the index of the second point
		 * @return boolean whether the points should be connected or no
		 */
    	bool differenceThreshold(sg::base::DataMatrix& testDataset, int i, int j);

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
    	sg::datadriven::LearnerTiming train(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes,
    						  const sg::base::RegularGridConfiguration& GridConfig, const sg::solver::SLESolverConfiguration& SolverConfig,
    						  const double lambda);

    	/**
    	 * TODO: Not implemented
    	 * @param testDataset
    	 * @param lambda
    	 * @return
    	 */
    	sg::datadriven::DMSystemMatrixBase* createDMSystem(sg::base::DataMatrix& testDataset, double lambda);

    	/**
    	 * Gets the last clustering result
    	 *
    	 * @return clustering result
    	 */
    	sg::base::DataVector getClusterAssignments();

    	/**
    	 * Precalculates the grid value for each point and stores the information a file.
    	 *
    	 * @param filename where to store
    	 * @param testDataset the data
    	 * @param GridConfig the grid configuration
    	 * @param SolverConfig the solver configuration
    	 * @param lamda the regularization parameter
    	 */
    	void precalculateGridValues(const char * filename, sg::base::DataMatrix& testDataset, const sg::base::RegularGridConfiguration& GridConfig, const sg::solver::SLESolverConfiguration& SolverConfig,
				  const double lamda);

    	/**
    	 * Load the precalculated values from a file
    	 * @param filename the file
    	 */
    	void loadPrecalculatedValues(const char * filename);

    	/**
    	 * Cluster function. The grid values must be calculated with ::precalculateGridValues and the cluster configuration must be set with the function ::setClusterConfiguration
    	 * @param trainDataset the data
    	 * @param GridConfig the grid configuration
    	 */
    	void cluster(sg::base::DataMatrix& trainDataset, const sg::base::RegularGridConfiguration& GridConfig);

    	/**
    	 * Set the cluster configuration
    	 *
    	 * @param ClusterConfig the cluster configuration
    	 */
    	void setClusterConfiguration(const sg::datadriven::DensityBasedClusteringConfiguration& ClusterConfig);

    	/**
    	 * Do a simple classification with the clustering result.The n biggest clusters form the classes.
    	 *
    	 * @param trainDataset the data set
    	 * @param components the clustering result
    	 * @param threshold
    	 * @return
    	 */
    	sg::base::DataVector postprocessing(sg::base::DataMatrix& trainDataset,sg::base::DataVector maxSize, int threshold);

    	/**
    	 * Prcalculates the n nearest neighbors and store the information in a file
    	 * @param filename the file
    	 * @param testDataset the data
    	 * @param n n nearest neighbors
    	 */
    	void precalculateNeighbors(const char * filename, sg::base::DataMatrix& testDataset, int n);

    	/**
    	 * Loads the precalculated neighbors
    	 * @param filename file
    	 */
    	void loadPrecalculatedNeighbors(const char * filename);
    };

  }

}

#endif /* LEARNERDENSITYCLUSTER_HPP */
