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

#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/Grid.hpp"

//#define ENABLE_DEBUG_MPI
#ifdef ENABLE_DEBUG_MPI
#define debugMPI(globalComm, messageStream) \
    { \
        int rank = globalComm->getMyRank(); \
        int rcvbuffer; \
        int token = 93453; \
        MPI_Barrier(MPI_COMM_WORLD); \
        if(rank > 0){ \
                MPI_Recv(&rcvbuffer, 1, MPI_INT, rank-1, token, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
        } \
        std::cout << "[" << rank << "] " <<  messageStream << std::endl; \
        std::cout.flush();\
        if(rank < globalComm->getNumRanks()-1){ \
            MPI_Send(&rank, 1, MPI_INT, rank+1, token, MPI_COMM_WORLD); \
        } \
        MPI_Barrier(MPI_COMM_WORLD); \
    }

#define debugMPI_0(messageStream) \
    { \
        int rank = globalComm->getMyRank(); \
        MPI_Barrier(MPI_COMM_WORLD); \
        if(rank == 0){ \
            std::cout << "[0] " <<  messageStream << std::endl; \
        } \
        MPI_Barrier(MPI_COMM_WORLD); \
    }
#else
#define debugMPI(globalComm, messageStream)
#define debugMPI_0(messageStream)
#endif

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

    void dataVectorAllToAll_alltoallv(sg::base::DataVector& alpha, int* distributionOffsets, int* distributionSizes);
    void dataVectorAllToAll_broadcasts(sg::base::DataVector& alpha, int* distributionOffsets, int* distributionSizes);
    void dataVectorAllToAll_sendreceive(sg::base::DataVector& alpha, int* distributionOffsets, int* distributionSizes);

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
     * broadcasts a specific range of a DataVector (depending on the rank of the process
     * and specified by distributionOffsets  and distributionSizes) to all other
     * processes.
     *
     * After this method has been executed, the entries of alpha from distributionOffsets[rank]
     * to distributionOffsets[rank] + distributionSizes[rank] - 1 of process with number rank
     * are available in all other processes in the same place. Overlapping regions lead to
     * undefined behaviour.
     *
     * example:
     *
     *      p1 p2 p3 p4                                    p1 p2 p3 p4
     * d[0] R  X  X  X                                d[0] R  R  R  R
     * d[1] S  X  X  X      distrOffsets: {0,4,2,5}   d[1] S  S  S  S
     * d[2] X  X  T  X     ------------------------>  d[2] T  T  T  T
     * d[3] X  X  U  X      distrSizes: {2,1,2,2}     d[3] U  U  U  U
     * d[4] X  V  X  X                                d[4] V  V  V  V
     * d[5] X  X  X  Y                                d[5] Y  Y  Y  Y
     * d[6] X  X  X  Z                                d[6] Z  Z  Z  Z
     *
     * @param alpha the DataVector to distribute
     * @param distributionOffsets array containing the offsets of data to distribute
     * @param distributionSizes array containing the sizes of data to distribute
     */
    void dataVectorAllToAll(sg::base::DataVector& alpha, int* distributionOffsets, int* distributionSizes);

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

    /**
     * calculate size and offset of fragment number rank for a domain with a size of totalSize for a distribution over procCount Fragments
     * distribute domain as equal as possible
     * (the difference between the minimum and maximum number of items is at most 1)
     * otherwise, a worker could have almost twice as much to do as all the others (example: total=127, proccount = 16)
     *
     * @param totalSize size of domain to distribute
     * @param procCount number of fragments (processes) to distribute domain into
     * @param rank specifies the number of the fragment for which to calculate size and offset
     * @param size output variable to put resulting size into
     * @param offset output variable to put resulting offset into
     */
    void calcDistributionFragment(int totalSize, int procCount, int rank, int *size, int *offset);

    /**
     * overloaded function.
     *
     * this function uses getNumRanks() and getMyRank() for the procCount and rank arguemnts
     *
     * @param totalSize size of domain to distribute
     * @param size output variable to put resulting size into
     * @param offset output variable to put resulting offset into
     */
    void calcDistributionFragment(int totalSize, int *size, int *offset);
};

}

}

#endif /* MPICOMMUNICATOR_HPP */
