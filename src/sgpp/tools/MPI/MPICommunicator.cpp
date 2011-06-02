/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "tools/MPI/MPICommunicator.hpp"

namespace sg
{
namespace parallel
{

MPICommunicator::MPICommunicator(int myid, int ranks) : myid_(myid), ranks_(ranks) { }

MPICommunicator::~MPICommunicator() { }

//void MPICommunicator::sendGridCoefficients(sg::base::DataVector& alpha, int dest_rank)
//{
//	MPI_Send((void*)alpha.getPointer(), (int)alpha.getSize(), MPI_DOUBLE, dest_rank, this->myid_, MPI_COMM_WORLD);
//}
//
//void MPICommunicator::broadcastGridCoefficients(sg::base::DataVector& alpha)
//{
//	for (int dest_rank = 1; dest_rank < this->ranks_; dest_rank++)
//	{
//		MPI_Send((void*)alpha.getPointer(), (int)alpha.getSize(), MPI_DOUBLE, dest_rank, this->myid_, MPI_COMM_WORLD);
//	}
//}

void MPICommunicator::broadcastGridCoefficientsFromRank0(sg::base::DataVector& alpha)
{
	MPI_Bcast((void*)alpha.getPointer(), (int)alpha.getSize(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

//void MPICommunicator::receiveGridCoefficients(sg::base::DataVector& alpha)
//{
//	MPI_Status status;
//
//	MPI_Recv((void*)alpha.getPointer(), (int)alpha.getSize(), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//}
//
//void MPICommunicator::aggregateGridCoefficients(sg::base::DataVector& alpha)
//{
//	for (int recv_rank = 1; recv_rank < this->ranks_; recv_rank++)
//	{
//		sg::base::DataVector tmp(alpha);
//		MPI_Status status;
//
//		MPI_Recv((void*)tmp.getPointer(), (int)tmp.getSize(), MPI_DOUBLE, recv_rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//
//		alpha.add(tmp);
//	}
//}

void MPICommunicator::reduceGridCoefficientsOnRank0(sg::base::DataVector& alpha)
{
	if (myid_ == 0)
	{
		MPI_Reduce(MPI_IN_PLACE, (void*)alpha.getPointer(), (int)alpha.getSize(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Reduce((void*)alpha.getPointer(), NULL, (int)alpha.getSize(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
}


void MPICommunicator::sendGrid(std::string& serialized_grid, int dest_rank)
{
	int nChars = (int)serialized_grid.length();
	MPI_Send(&nChars, 1, MPI_INT, dest_rank, this->myid_, MPI_COMM_WORLD);
	MPI_Send((void*)serialized_grid.c_str(), (int)serialized_grid.length(), MPI_CHAR, dest_rank, this->myid_, MPI_COMM_WORLD);
}

void MPICommunicator::broadcastGrid(std::string& serialized_grid)
{
	int nChars = (int)serialized_grid.length();
	for (int dest_rank = 1; dest_rank < this->ranks_; dest_rank++)
	{
		MPI_Send(&nChars, 1, MPI_INT, dest_rank, this->myid_, MPI_COMM_WORLD);
		MPI_Send((void*)serialized_grid.c_str(), (int)serialized_grid.length(), MPI_CHAR, dest_rank, this->myid_, MPI_COMM_WORLD);
	}
}

void MPICommunicator::receiveGrid(std::string& serialized_grid)
{
	char* recv_buffer;
	int count;
	MPI_Status status;

	MPI_Recv((void*)(&count), 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	recv_buffer = new char[count];
	MPI_Recv((void*)recv_buffer, count, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

	serialized_grid = "";
	serialized_grid.assign(recv_buffer);

	delete recv_buffer;
}

void MPICommunicator::broadcastGridStorage(std::string& serialized_grid_storage)
{
	int nChars = (int)serialized_grid_storage.length();
	for (int dest_rank = 1; dest_rank < this->ranks_; dest_rank++)
	{
		MPI_Send(&nChars, 1, MPI_INT, dest_rank, this->myid_, MPI_COMM_WORLD);
		MPI_Send((void*)serialized_grid_storage.c_str(), (int)serialized_grid_storage.length(), MPI_CHAR, dest_rank, this->myid_, MPI_COMM_WORLD);
	}
}

void MPICommunicator::receiveGridStorage(std::string& serialized_grid_storage)
{
	char* recv_buffer;
	int count;
	MPI_Status status;

	MPI_Recv((void*)(&count), 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	recv_buffer = new char[count];
	MPI_Recv((void*)recv_buffer, count, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

	serialized_grid_storage = "";
	serialized_grid_storage.assign(recv_buffer);

	delete recv_buffer;
}

void MPICommunicator::Barrier()
{
	MPI_Barrier(MPI_COMM_WORLD);
}

void MPICommunicator::Abort()
{
	MPI_Abort(MPI_COMM_WORLD, -1);
}

//void MPICommunicator::broadcastControl(char ctrl)
//{
//	for (int dest_rank = 1; dest_rank < this->ranks_; dest_rank++)
//	{
//		MPI_Send((void*)&ctrl, 1, MPI_CHAR, dest_rank, this->myid_, MPI_COMM_WORLD);
//	}
//}
//
//char MPICommunicator::receiveControl()
//{
//	MPI_Status status;
//	char result;
//
//	MPI_Recv((void*)&result, 1, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//
//	return result;
//}

void MPICommunicator::broadcastControlFromRank0(char* ctrl)
{
	MPI_Bcast((void*)ctrl, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
}

int MPICommunicator::getMyRank()
{
	return this->myid_;
}

int MPICommunicator::getNumRanks()
{
	return this->ranks_;
}

}
}
