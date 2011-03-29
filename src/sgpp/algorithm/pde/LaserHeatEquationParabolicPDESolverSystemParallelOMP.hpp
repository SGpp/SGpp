/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef LASERHEATEQUATIONPARABOLICPDESOLVERSYSTEMPARALLELOMP_HPP
#define LASERHEATEQUATIONPARABOLICPDESOLVERSYSTEMPARALLELOMP_HPP

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
class LaserHeatEquationParabolicPDESolverSystemParallelOMP : public HeatEquationParabolicPDESolverSystemParallelOMP
{
private:
	/// the velocity of the rotating laser beam
	double beam_velocity_;
	/// the laser beam's expansion
	double heat_sigma_;
	/// the max. refinement level
	size_t max_level_;

	/// current offset of laser position (x direction)
	double laser_x_offset_;
	/// position of laser at beginning (x direction)
	double laser_x_start_;
	/// position of laser in last time step (x direction)
	double laser_x_last_;
	/// indicator of sign for y offset
	size_t y_sign_switch;

public:
	/**
	 * Std-Constructor
	 *
	 * @param beam_velocity the velocity of the rotating laser beam
	 * @param heat_sigma the laser beam's expansion
	 * @param max_level the max. refinement level
	 * @param SparseGrid reference to the sparse grid
	 * @param alpha the sparse grid's coefficients
	 * @param a the heat coefficient
	 * @param TimestepSize the size of one timestep used in the ODE Solver
	 * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for explicit Euler,
	 *  							ImEul for implicit Euler, CrNic for Crank Nicolson solver
	 */
	LaserHeatEquationParabolicPDESolverSystemParallelOMP(double beam_velocity, double heat_sigma, size_t max_level, Grid& SparseGrid, DataVector& alpha, double a, double TimestepSize, std::string OperationMode = "ExEul");

	/**
	 * Std-Destructor
	 */
	virtual ~LaserHeatEquationParabolicPDESolverSystemParallelOMP();

	void finishTimestep(bool isLastTimestep = false);

	void startTimestep();
};

}

#endif /* LASERHEATEQUATIONPARABOLICPDESOLVERSYSTEMPARALLELOMP_HPP */
