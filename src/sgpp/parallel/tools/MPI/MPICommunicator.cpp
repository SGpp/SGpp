/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

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
	int numRanks = getNumRanks();
	int myRank = getMyRank();
	int mySendSize = distributionSizes[myRank];
	int mySendOffset = distributionOffsets[myRank];

	int *sendSizes = new int[numRanks];
	int *sendOffsets = new int[numRanks];
	for(int i = 0; i<numRanks;i++){
		sendSizes[i] = mySendSize;
		sendOffsets[i] = mySendOffset;
	}

	sg::base::DataVector tmp(alpha.getSize());
	MPI_Alltoallv(alpha.getPointer(), sendSizes, sendOffsets, MPI_DOUBLE,
				  tmp.getPointer(), distributionSizes, distributionOffsets, MPI_DOUBLE, MPI_COMM_WORLD);
	alpha.copyFrom(tmp);

	delete[] sendSizes;
	delete[] sendOffsets;
}

void MPICommunicator::dataVectorAllToAll(base::DataVectorSP &alpha, int *distributionOffsets, int *distributionSizes)
{
	int numRanks = getNumRanks();
	int myRank = getMyRank();
	int mySendSize = distributionSizes[myRank];
	int mySendOffset = distributionOffsets[myRank];

	int *sendSizes = new int[numRanks];
	int *sendOffsets = new int[numRanks];
	for(int i = 0; i<numRanks;i++){
		sendSizes[i] = mySendSize;
		sendOffsets[i] = mySendOffset;
	}

	sg::base::DataVectorSP tmp(alpha.getSize());
	MPI_Alltoallv(alpha.getPointer(), sendSizes, sendOffsets, MPI_FLOAT,
				  tmp.getPointer(), distributionSizes, distributionOffsets, MPI_FLOAT, MPI_COMM_WORLD);
	alpha.copyFrom(tmp);

	delete[] sendSizes;
	delete[] sendOffsets;
}

void MPICommunicator::IsendToAll(double *ptr, size_t size, int tag)
{
	MPI_Request* req = new MPI_Request();
	for(int rank = 0; rank<getNumRanks(); rank++){
		if (rank == getMyRank()){
			continue;
		}

		MPI_Isend(ptr, size, MPI_DOUBLE, rank, tag, MPI_COMM_WORLD, req);
	}
	delete req;
}

void MPICommunicator::IrecvFromAll(double *ptr, int *global_sizes, int *global_offsets, int *sizes, int *offsets, int *tag, int chunkCount, MPI_Request *dataRecvRequests)
{
	int rank;
	//posting temp reveices
	rank = 0;
	for(int i = 0; i<chunkCount; i++){
		// adjust rank to match chunk i
		if(i>=global_offsets[rank] + global_sizes[rank]){
			rank++;
		}
		// skip segments of this process, they are already there
		if(rank == getMyRank()){
			//i=global_offsets[rank] + global_sizes[rank]-1;
			dataRecvRequests[i] = MPI_REQUEST_NULL;
			// continue does the i++ (like after every iteration), so we are
			// at _mpi_data_offsets_global[rank] + _mpi_data_sizes_global[rank] at
			// the beginning of the next iteration, which means that we skipped rank mpi_myrank
			continue;
		}
		MPI_Irecv(&ptr[offsets[i]], sizes[i], MPI_DOUBLE, rank, tag[i], MPI_COMM_WORLD, &dataRecvRequests[i]);
	}
}

void MPICommunicator::putToAll(double* ptr, int winOffset, int count, MPI_Win win)
{
	for(int i = 0; i<getNumRanks(); i++){
		if(i==getMyRank()){
			//continue;
		}
		MPI_Put(ptr, count, MPI_DOUBLE, i, winOffset, count, MPI_DOUBLE, win);
	}
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
    int numRanks = getNumRanks();
    int myRank = getMyRank();
    int mySendSize = distributionSizes[myRank];
    int mySendOffset = distributionOffsets[myRank];

    debugMPI(this, "read distribution vars");

    int *sendSizes = new int[numRanks];
    int *sendOffsets = new int[numRanks];
    for(int i = 0; i<numRanks;i++){
        sendSizes[i] = mySendSize;
        sendOffsets[i] = mySendOffset;
    }

    debugMPI(this, "initialized sendVars");

//    std::ostringstream strstr;

//    strstr << "sendSizes: ";
//    for(int i = 0; i<numRanks; i++){
//        strstr << sendSizes[i] << " - ";
//    }
//    debugMPI(this, strstr.str());
//    strstr.str("");

//    strstr << "sendOffsets: ";
//    for(int i = 0; i<numRanks; i++){
//        strstr << sendOffsets[i] << " - ";
//    }
//    debugMPI(this, strstr.str());
//    strstr.str("");

//    strstr << "distributionSizes: ";
//    for(int i = 0; i<numRanks; i++){
//        strstr << distributionSizes[i] << " - ";
//    }
//    debugMPI(this, strstr.str());
//    strstr.str("");

//    strstr << "distributionOffsets: ";
//    for(int i = 0; i<numRanks; i++){
//        strstr << distributionOffsets[i] << " - ";
//    }
//    debugMPI(this, strstr.str());
//    strstr.str("");

//    debugMPI(this, strstr.str());
//    strstr.str("");


    sg::base::DataVector tmp(alpha.getSize());
    MPI_Alltoallv(alpha.getPointer(), sendSizes, sendOffsets, MPI_DOUBLE,
                  tmp.getPointer(), distributionSizes, distributionOffsets, MPI_DOUBLE, MPI_COMM_WORLD);
    alpha.copyFrom(tmp);


    // this seems to work, altough the standard does not guarantee this,
    // also, the standard is inconsistent regarding the use of MPI_IN_PLACE with MPI_Alltoallv()
//    MPI_Alltoallv(alpha.getPointer(), sendSizes, sendOffsets, MPI_DOUBLE,
//                  alpha.getPointer(), distributionSizes, distributionOffsets, MPI_DOUBLE, MPI_COMM_WORLD);

    //the following does not work, MPI_IN_PLACE is not allowed
//    MPI_Alltoallv(MPI_IN_PLACE, distributionSizes, distributionSizes, MPI_DOUBLE,
//                  alpha.getPointer(), distributionSizes, distributionOffsets, MPI_DOUBLE, MPI_COMM_WORLD);

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
