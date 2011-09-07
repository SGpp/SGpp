/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef MPICOMMUNICATOR_HPP
#define MPICOMMUNICATOR_HPP

#include <mpi.h>

#include <string>

#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"
#include "grid/Grid.hpp"

namespace sg
{
namespace parallel
{

/**
 * This class provides standard tasks, like
 * sending and receiving coefficients and grids,
 * needed when dealing with distributed
 * sparse grid applications.
 *
 * @version $HEAD$
 */
class MPICommunicator
{
private:
	/// The rank that owns this instance
	int myid_;
	/// Number of ranks in the whole system
	int ranks_;

public:
	/*
	 * std-constructor
	 *
	 * @param myid current rank's id
	 * @param ranks number of ranks in the team
	 */
	MPICommunicator(int myid, int ranks);

	/*
	 * std-destructor
	 */
	~MPICommunicator();

	/**
	 * sends grid coefficients to all ranks greater zero
	 * here MPI_Bcast is used.
	 *
	 * @param alpha grid coefficients that should be sent
	 */
	void broadcastGridCoefficientsFromRank0(sg::base::DataVector& alpha);

	/**
	 * Reduces the grid coefficients on rank 0 using
	 * MPI's reduce routine
	 *
	 * @param alpha sg::base::DataVector to which all other rank's grid coefficients should be added
	 */
	void reduceGridCoefficientsOnRank0(sg::base::DataVector& alpha);

	/**
	 * Reduces the grid coefficients on of all ranks and all ranks using
	 * MPI's reduce routine
	 *
	 * @param alpha sg::base::DataVector to which all other rank's grid coefficients should be added
	 */
	void reduceGridCoefficients(sg::base::DataVector& alpha);

	/**
	 * sends a serialized grid to a specific rank
	 *
	 * @param serialized_grid serialized grid that should be sent
	 * @param dest_rank destination rank
	 */
	void sendGrid(std::string& serialized_grid, int dest_rank);

	/**
	 * sends a serialized grid to all ranks
	 *
	 * @param serialized_grid serialized grid that should be sent
	 */
	void broadcastGrid(std::string& serialized_grid);

	/**
	 * receives a serialized grid an stores it into serialized_grid
	 *
	 * @param serialized_grid the serialized grid
	 */
	void receiveGrid(std::string& serialized_grid);

	/**
	 * sends a serialized grid storage to all ranks
	 *
	 * @param serialized_grid_storage serialized grid storage that should be sent
	 */
	void broadcastGridStorage(std::string& serialized_grid_storage);

	/**
	 * receives a serialized grid storage an stores it into serialized_grid_storage
	 *
	 * @param serialized_grid_storage the serialized grid storage
	 */
	void receiveGridStorage(std::string& serialized_grid_storage);

	/**
	 * sends a control character to all ranks greater the zero
	 * MPI_Bcast is used.
	 *
	 * @param ctrl Control character that should be sent
	 */
	void broadcastControlFromRank0(char* ctrl);

	/**
	 * Implements a Barrier for all tasks
	 */
	void Barrier();

	/**
	 * Implements a Barrier for all tasks
	 */
	void Abort();

	/**
	 * @return return the rank id of the own of this class
	 */
	int getMyRank();

	/**
	 * @return returns the number of MPI tasks in the parallel environment
	 */
	int getNumRanks();
};

}

}

#endif /* MPICOMMUNICATOR_HPP */
