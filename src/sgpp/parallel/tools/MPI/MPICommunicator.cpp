/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/tools/MPI/MPICommunicator.hpp"

namespace sg
{
namespace parallel
{

MPICommunicator::MPICommunicator(int myid, int ranks) : myid_(myid), ranks_(ranks) { }

MPICommunicator::~MPICommunicator() { }

void MPICommunicator::broadcastGridCoefficientsFromRank0(sg::base::DataVector& alpha)
{
    MPI_Bcast((void*)alpha.getPointer(), (int)alpha.getSize(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

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

void MPICommunicator::reduceGridCoefficients(sg::base::DataVector& alpha)
{
	MPI_Allreduce(MPI_IN_PLACE, (void*)alpha.getPointer(), (int)alpha.getSize(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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

void MPICommunicator::broadcastControlFromRank0(char* ctrl)
{
    MPI_Bcast((void*)ctrl, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
}

void MPICommunicator::dataVectorAllToAll(base::DataVector &alpha, int *distributionOffsets, int *distributionSizes)
{
    dataVectorAllToAll_alltoallv(alpha, distributionOffsets, distributionSizes);
}

int MPICommunicator::getMyRank()
{
	return this->myid_;
}

int MPICommunicator::getNumRanks()
{
    return this->ranks_;
}

void MPICommunicator::dataVectorAllToAll_alltoallv(base::DataVector &alpha, int *distributionOffsets, int *distributionSizes)
{
    debugMPI(this, "starting alltoall");
    int numRanks = getNumRanks();
    int myRank = getMyRank();
    int mySendSize = distributionSizes[myRank];
    int mySendOffset = distributionOffsets[myRank];

    double *alphaPointer = alpha.getPointer();

    debugMPI(this, "read distribution vars");

    int *sendSizes = new int[numRanks];
    int *sendOffsets = new int[numRanks];
    for(int i = 0; i<numRanks;i++){
        sendSizes[i] = mySendSize;
        sendOffsets[i] = mySendOffset;
    }

    debugMPI(this, "initialized sendVars");

    sg::base::DataVector tmp(alpha.getSize());
    //tmp = new sg::base::DataVector(alpha.getSize());
    //sg::base::DataVector tmp(alpha);

    debugMPI(this, "before communication - alltoall; vector size: " << alpha.getSize());

    std::ostringstream strstr;

    strstr << "sendSizes: ";
    for(int i = 0; i<numRanks; i++){
        strstr << sendSizes[i] << " - ";
    }
    debugMPI(this, strstr.str());
    strstr.str("");

    strstr << "sendOffsets: ";
    for(int i = 0; i<numRanks; i++){
        strstr << sendOffsets[i] << " - ";
    }
    debugMPI(this, strstr.str());
    strstr.str("");

    strstr << "distributionSizes: ";
    for(int i = 0; i<numRanks; i++){
        strstr << distributionSizes[i] << " - ";
    }
    debugMPI(this, strstr.str());
    strstr.str("");

    strstr << "distributionOffsets: ";
    for(int i = 0; i<numRanks; i++){
        strstr << distributionOffsets[i] << " - ";
    }
    debugMPI(this, strstr.str());
    strstr.str("");

    strstr << "some values: ";
    for(int i = 0; i<10; i++){
        strstr << alphaPointer[i] << "; ";
    }
    for(int i = alpha.getSize()-10; i<alpha.getSize(); i++){
        strstr << alphaPointer[i] << "; ";
    }
    debugMPI(this, strstr.str());
    strstr.str("");


    MPI_Alltoallv(alpha.getPointer(), sendSizes, sendOffsets, MPI_DOUBLE,
                  tmp.getPointer(), distributionSizes, distributionOffsets, MPI_DOUBLE, MPI_COMM_WORLD);
    alpha.copyFrom(tmp);


//    MPI_Alltoallv(alphaPointer, sendSizes, sendOffsets, MPI_DOUBLE,
//                  alphaPointer, distributionSizes, distributionOffsets, MPI_DOUBLE, MPI_COMM_WORLD);

//    MPI_Alltoallv(MPI_IN_PLACE, distributionSizes, distributionSizes, MPI_DOUBLE,
//                  alphaPointer, distributionSizes, distributionOffsets, MPI_DOUBLE, MPI_COMM_WORLD);

    debugMPI(this, "after communication - alltoallv; ");

    delete[] sendSizes;
    delete[] sendOffsets;
}

void MPICommunicator::dataVectorAllToAll_broadcasts(base::DataVector &alpha, int *distributionOffsets, int *distributionSizes)
{
    double *data = alpha.getPointer();
    for(int i = 0; i<getNumRanks(); i++){
        MPI_Bcast(&data[distributionOffsets[i]], distributionSizes[i], MPI_DOUBLE, i, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void MPICommunicator::dataVectorAllToAll_sendreceive(base::DataVector &alpha, int *distributionOffsets, int *distributionSizes)
{
    double *data = alpha.getPointer();



    MPI_Request reqs[getNumRanks()];
    int myRank = getMyRank();
    for(int otherProc = 0; otherProc < getNumRanks(); ++otherProc){
        if (otherProc == myRank) {
            continue;
        }
        MPI_Request r;
        MPI_Isend(&data[distributionOffsets[myRank]], distributionSizes[myRank], MPI_DOUBLE, otherProc, 0, MPI_COMM_WORLD, &r);
        MPI_Request_free(&r);
    }

    std::cout << "["<< myRank <<"] after isend; " << std::endl;


    for(int otherProc = 0; otherProc < getNumRanks(); ++otherProc) {
        std::cout << "["<< myRank <<"] before irecv; otherp:" << otherProc << std::endl;
        if (otherProc == myRank) {
            reqs[otherProc] = MPI_REQUEST_NULL;
            continue;
        }
        if(myRank == 2) {
            std::cout << "["<< myRank <<"]  offset: " << distributionOffsets[otherProc] << " sizes:" << distributionSizes[otherProc] << std::endl;
        }
        MPI_Irecv(&data[distributionOffsets[otherProc]], distributionSizes[otherProc], MPI_DOUBLE, otherProc, 0, MPI_COMM_WORLD, &reqs[otherProc]);
        std::cout << "["<< myRank <<"] after irecv; otherp:" << otherProc << std::endl;
    }

    MPI_Waitall(getNumRanks(), reqs, MPI_STATUSES_IGNORE);
    std::cout << "["<< myRank <<"] after waitall; " << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
}

}

}
