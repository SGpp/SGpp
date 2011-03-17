/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef PDESOLVER_HPP
#define PDESOLVER_HPP

#include "grid/Grid.hpp"
#include "grid/common/BoundingBox.hpp"

#include "data/DataVector.hpp"

#include "tools/finance/IOToolBonnSG.hpp"
#include "tools/common/GridPrinter.hpp"

#include <vector>

namespace sg
{

/**
 * This class provides defines a implements basic task and tools which are
 * useful while solving PDEs. E.g. grid construction, simple grid evaluation tools
 * grid printing, initial grid refinement etc.
 *
 * @version $HEAD$
 */
class PDESolver
{
protected:
	/// The Sparse Grid needed in this classificator
	Grid* myGrid;
	/// the number of levels used for an regular grid
	size_t levels;
	/// the dimension of the grid
	size_t dim;
	/// stores if the grid was created inside the solver
	bool bGridConstructed;
	/// Stores Pointer to the Grid's Bounding Box
	BoundingBox* myBoundingBox;
	/// Stores Pointer to the Girs's Storage
	GridStorage* myGridStorage;

	/**
	 * This function calculates for every grid point the value
	 * of a normal distribution given by norm_mu and norm_sigma.
	 * The result is stored dehierarchized in alpha.
	 *
	 * @param alpha contains dehierarchized sparse grid coefficients containing the values of the multi dimensional normal distribution after call
	 * @param norm_mu the expected values of the normal distribution for every grid dimension
	 * @param norm_sigma the standard deviation of the normal distribution for every grid dimension
	 */
	virtual void getGridNormalDistribution(DataVector& alpha, std::vector<double>& norm_mu, std::vector<double>& norm_sigma);

public:
	/**
	 * Std-Constructor of the solver
	 */
	PDESolver();

	/**
	 * Std-Destructor of the solver
	 */
	virtual ~PDESolver();

	/**
	 * Use this routine the construct a regular grid to solve a PDE
	 *
	 * @param myBoundingBox reference to a bounding box that describes the grid
	 * @param level number of the regular's grid levels
	 */
	virtual void constructGrid(BoundingBox& myBoundingBox, size_t level) = 0;

	/**
	 * Sets the grid used in this BlackScholes Solver by an given serialized string
	 * of the grid.
	 *
	 * @param serializedGrid a string that describes the grid that should be used in this solver
	 */
	void setGrid(const std::string& serializedGrid);

	/**
	 * gets the a string the describes the grid which is currently used to solve
	 *
	 * @return a string containing a serialized grid
	 */
	std::string getGrid() const;

	/**
	 * deletes the grid created within that solver
	 */
	void deleteGrid();

	/**
	 * Refines a grid by taking the grid's coefficients into account. This refinement method
	 * refines the grid based on the surplus by refining grid points with big surpluses
	 * first. The number of grid points to refine may be specified by the numRefinePoints parameter.
	 *
	 * @param alpha a DataVector containing the grids coefficients
	 * @param numRefinePoints the number of grid points that should be refined; if this smaller than zero -> all refineable points will be refined
	 * @param dThreshold Threshold for a point's surplus for refining this point
	 */
	void refineInitialGridSurplus(DataVector& alpha, int numRefinePoints, double dThreshold);

	/**
	 * Refines a grid by taking the grid's coefficients into account. This refinement method
	 * refines the grid based on the surplus by refining grid points with big surpluses
	 * first.
	 * The grid is refined to max. Level!
	 *
	 * @param alpha a DataVector containing the grids coefficients
	 * @param dThreshold Threshold for a point's surplus for refining this point
	 * @param maxLevel maxLevel of refinement
	 */
	void refineInitialGridSurplusToMaxLevel(DataVector& alpha, double dThreshold, size_t maxLevel);

	/**
	 * Refines a grid by taking the grid's coefficients into account. This refinement method
	 * refines the grid based on the surplus by refining grid points with big surpluses
	 * first. The number of grid points to refine may be specified by the numRefinePoints parameter.
	 *
	 * This functions refines the grid only in subdomain specified by a multi-dimensional
	 * normal distribution. The normal distribution is given by the parameters norm_mu
	 * and norm_sigma which are d-dimensional vectors.
	 *
	 * @param alpha a DataVector containing the grids coefficients
	 * @param numRefinePoints the number of grid points that should be refined; if this smaller than zero -> all refineable points will be refined
	 * @param dThreshold Threshold for a point's surplus for refining this point
	 * @param norm_mu the expected values of the normal distribution for every grid dimension
	 * @param norm_sigma the standard deviation of the normal distribution for every grid dimension
	 */
	void refineInitialGridSurplusSubDomain(DataVector& alpha, int numRefinePoints, double dThreshold, std::vector<double>& norm_mu, std::vector<double>& norm_sigma);

	/**
	 * Refines a grid by taking the grid's coefficients into account. This refinement method
	 * refines the grid based on the surplus by refining grid points with big surpluses
	 * first.
	 * The grid is refined to max. Level!
	 *
	 * This functions refines the grid only in subdomain specified by a multi-dimensional
	 * normal distribution. The normal distribution is given by the parameters norm_mu
	 * and norm_sigma which are d-dimensional vectors.
	 *
	 * @param alpha a DataVector containing the grids coefficients
	 * @param dThreshold Threshold for a point's surplus for refining this point
	 * @param maxLevel maxLevel of refinement
	 * @param norm_mu the expected values of the normal distribution for every grid dimension
	 * @param norm_sigma the standard deviation of the normal distribution for every grid dimension
	 */
	void refineInitialGridSurplusToMaxLevelSubDomain(DataVector& alpha, double dThreshold, size_t maxLevel, std::vector<double>& norm_mu, std::vector<double>& norm_sigma);

	/**
	 * Coarsens a grid by taking the grid's coefficients into account. This coarsen method
	 * coarsens the grid based on the surplus by coarsening grid points with small surpluses
	 * first.
	 *
	 * @param alpha a DataVector containing the grids coefficients
	 * @param dThreshold Threshold for a point's surplus for coarsening this point
	 */
	void coarsenInitialGridSurplus(DataVector& alpha, double dThreshold);

	/**
	 * Use this routine the construct a regular grid to solve the multi-dimensional Black Scholes Equation
	 *
	 * Use this routine if you want to solve a problem stored in the format provided by the solving system
	 * released by the University of Bonn, Germany
	 *
	 * @param tfilename absolute path of the file
	 * @param emptyAlpha reference to a DataVector object that contains no elements
	 * @param ishierarchized is set to true if alpha contains surplus after reading the file, otherwise false
	 */
	void constructGridBonn(std::string tfilename, DataVector& emptyAlpha, bool& ishierarchized);

	/**
	 * Use this routine if you wnat to store a grid in the format provided by the solving system
	 * released by the University of Bonn, Germany
	 *
	 * @param tfilename absolute path of the file
	 * @param alpha reference to a DataVector object that contains the gird ansatzfunction's coefficients
	 * @param ishierarchized set to true, the export is done on the nodal basis
	 */
	void storeGridBonn(std::string tfilename, DataVector& alpha, bool ishierarchized);

	/**
	 * Determines the value of the function in the d-dimensional space
	 *
	 * @param evalPoint coordinates of the point at which the function should be evaluated
	 * @param alpha the ansatzfunctions' coefficients
	 *
	 * @return price of option for given point
	 */
	double evaluatePoint(std::vector<double>& evalPoint, DataVector& alpha);

	/**
	 * Evaluates the sparse grid's function given by the stored grid and the alpha coefficients.
	 * on different points specified in EvaluationPoints and stores the result into FunctionValues.
	 *
	 * @param alpha the sparse grid's coefficients
	 * @param FunctionValues DataVector into the which the result of function's evaluation is stored
	 * @param EvaluationPoints DataMatrix that contains the points at which the sparse grid's function is evaluated
	 */
	void evaluateCuboid(DataVector& alpha, DataVector& FunctionValues, DataMatrix& EvaluationPoints);

	/**
	 * This is some kind of debug functionality. It writes a file,
	 * that can be used with gnuplot the print the grid.
	 *
	 * Is only implemented for 1D and 2D grids!
	 *
	 * @param alpha the coefficients of the Sparse Gird's basis functions
	 * @param PointesPerDimension the distance between evaluation points
	 * @param tfilename absolute path to file into which the grid's evaluation is written
	 */
	void printGrid(DataVector& alpha, double PointesPerDimension, std::string tfilename) const;

	/**
	 * This is some kind of debug functionality. It writes a file,
	 * that can be used with gnuplot the print the grid.
	 *
	 * Is only implemented for 2D grids!
	 *
	 * @param alpha the coefficients of the Sparse Gird's basis functions
	 * @param PointesPerDimension the distance between evaluation points
	 * @param GridArea the area in which the function should be plotted
	 * @param tfilename absolute path to file into which the grid's evaluation is written
	 */
	void printGridDomain(DataVector& alpha, double PointesPerDimension, BoundingBox& GridArea, std::string tfilename) const;

	/**
	 * Prints the Grid Points of the Sparse Grid either with their node basis value
	 * or their hierarchical surplus
	 *
	 * This function is available for all dimensions
	 *
	 * @param alpha the coefficients of the grid's ansatzfunctions
	 * @param tfilename absoulte path to the file the grid is written into
	 * @param bSurplus specifies whether the surplus (true) or the node basis value (false) is written
	 */
	void printSparseGrid(DataVector& alpha, std::string tfilename, bool bSurplus) const;

	/**
	 * Prints the Grid Points of the Sparse Grid either with their node basis value
	 * or their hierarchical surplus
	 *
	 * This function is available for all dimensions.
	 *
	 * The coordinates of the grid points are pushed the exp function. So
	 * log transformed grids can be plotted in cartesion coordinates.
	 *
	 * @param alpha the coefficients of the grid's ansatzfunctions
	 * @param tfilename absoulte path to the file the grid is written into
	 * @param bSurplus specifies whether the surplus (true) or the node basis value (false) is written
	 */
	void printSparseGridExpTransform(DataVector& alpha, std::string tfilename, bool bSurplus) const;

	/**
	 * use this to determine the number of grid points, used to solve
	 * the current problem
	 *
	 * @return the number of grid points
	 */
	size_t getNumberGridPoints() const;

	/**
	 * use this to determine the number of inner grid points, used to solve
	 * the current problem
	 *
	 * @return the number of inner grid points
	 */
	size_t getNumberInnerGridPoints() const;

	/**
	 * use this the determine the number of dimensions that are currently used
	 * in the solver.
	 *
	 * @return returns the number of the grid's dimensions, if the grid isn't constructed, yet it returns 0
	 */
	size_t getNumberDimensions() const;
};

}

#endif /* PDESOLVER_HPP */
