/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef LASERHEATEQUATIONPARABOLICPDESOLVERSYSTEMPARALLELOMP2D_HPP
#define LASERHEATEQUATIONPARABOLICPDESOLVERSYSTEMPARALLELOMP2D_HPP

#include "algorithm/pde/HeatEquationParabolicPDESolverSystemParallelOMP.hpp"

namespace sg
{

/**
 * This class implements the ParabolicPDESolverSystem for the
 * Heat Equation.
 *
 * In this case the Heat Equation is used to solve a rotating laser
 * beam.
 */
class LaserHeatEquationParabolicPDESolverSystemParallelOMP2D : public HeatEquationParabolicPDESolverSystemParallelOMP
{
private:
	/// the velocity of the rotating laser beam
	double beam_velocity_;
	/// the laser beam's expansion
	double heat_sigma_;
	/// the max. refinement level
	size_t max_level_;
	/// heating the grid initialization
	double heat_;
	/// threshold for refinening during solution
	double refine_threshold_;
	/// threshold for coarsening during solution
	double coarsen_threshold_;
	/// number of calculated timesteps
	size_t done_steps_;

public:
	/**
	 * Std-Constructor
	 *
	 * @param beam_velocity the velocity of the rotating laser beam
	 * @param heat_sigma the laser beam's expansion
	 * @param max_level the max. refinement level
	 * @param heat heating the grid initialization
	 * @param refine_threshold threshold for refinening during solution
	 * @param coarsen_threshold threshold for coarsening during solution
	 * @param SparseGrid reference to the sparse grid
	 * @param alpha the sparse grid's coefficients
	 * @param a the heat coefficient
	 * @param TimestepSize the size of one timestep used in the ODE Solver
	 * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for explicit Euler,
	 *  							ImEul for implicit Euler, CrNic for Crank Nicolson solver
	 */
	LaserHeatEquationParabolicPDESolverSystemParallelOMP2D(double beam_velocity, double heat_sigma, size_t max_level, double heat, double refine_threshold, double coarsen_threshold, Grid& SparseGrid, DataVector& alpha, double a, double TimestepSize, std::string OperationMode = "ExEul");

	/**
	 * Std-Destructor
	 */
	virtual ~LaserHeatEquationParabolicPDESolverSystemParallelOMP2D();

	void finishTimestep(bool isLastTimestep = false);

	void startTimestep();
};

}

#endif /* LASERHEATEQUATIONPARABOLICPDESOLVERSYSTEMPARALLELOMP2D_HPP */
