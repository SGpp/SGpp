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
	double laser_velo_;

public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param alpha the sparse grid's coefficients
	 * @param a the heat coefficient
	 * @param TimestepSize the size of one timestep used in the ODE Solver
	 * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for explicit Euler,
	 *  							ImEul for implicit Euler, CrNic for Crank Nicolson solver
	 */
	LaserHeatEquationParabolicPDESolverSystemParallelOMP(Grid& SparseGrid, DataVector& alpha, double a, double TimestepSize, std::string OperationMode = "ExEul");

	/**
	 * Std-Destructor
	 */
	virtual ~LaserHeatEquationParabolicPDESolverSystemParallelOMP();

	void finishTimestep(bool isLastTimestep = false);

	void startTimestep();
};

}

#endif /* LASERHEATEQUATIONPARABOLICPDESOLVERSYSTEMPARALLELOMP_HPP */
