/*
 * Task.cpp
 *
 *  Created on: Jul 10, 2014
 *      Author: heenemo
 */

#include "../../distributedcombigrid/task/Task.hpp"

#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace combigrid {

Task::Task() {
}

Task::Task(DimType dim, LevelVector& l, LevelVector& nmax, LevelVector& lmin,
    std::vector<bool>& boundary, real coeff, LoadModel* loadModel) :
    dim_(dim), l_(l), boundary_(boundary), id_(count++), loadModel_(
        loadModel), isFinished_(false) {
  assert(dim_ > 0);
  assert(l_.size() == dim_);
}

Task::~Task() {
}

int Task::count = 0;

void Task::send(Task** t, RankType dst, CommunicatorType& comm) {
  // save data to archive
  std::stringstream ss;
  {
    boost::archive::text_oarchive oa(ss);
    // write class instance to archive
    oa << *t;
  }
  // create mpi buffer of archive
  std::string s = ss.str();
  int bsize = static_cast<int>(s.size());
  char* buf = const_cast<char*>(s.c_str());
  MPI_Send(buf, bsize, MPI_CHAR, dst, 0, comm);
}

void Task::receive(Task** t, RankType src, CommunicatorType& comm) {
  // receive size of message
  // todo: not really necessary since size known at compile time
  MPI_Status status;
  int bsize;
  MPI_Probe(src, 0, comm, &status);
  MPI_Get_count(&status, MPI_CHAR, &bsize);

  // create buffer of appropriate size and receive
  std::vector<char> buf(bsize);
  MPI_Recv(&buf[0], bsize, MPI_CHAR, src, 0, comm, MPI_STATUS_IGNORE);

  // create and open an archive for input
  std::string s(&buf[0], bsize);
  std::stringstream ss(s);
  {
    boost::archive::text_iarchive ia(ss);
    // read class state from archive
    ia >> *t;
  }
}

void Task::broadcast(Task** t, RankType root, CommunicatorType& comm) {
  RankType myID;
  MPI_Comm_rank(comm, &myID);

  char* buf = NULL;
  int bsize;

  // root writes object data into buffer
  std::string s;
  if (myID == root) {
    // save data to archive
    std::stringstream ss;
    {
      boost::archive::text_oarchive oa(ss);
      // write class instance to archive
      oa << *t;
    }
    // create mpi buffer of archive
    s = ss.str();
    bsize = static_cast<int>(s.size());
    buf = const_cast<char*>(s.c_str());
  }

  // root broadcasts object size
  MPI_Bcast(&bsize, 1, MPI_INT, root, comm);

  // non-root procs create buffer which is large enough
  std::vector<char> tmp(bsize);
  if (myID != root) {
    buf = &tmp[0];
  }

  // broadcast of buffer
  MPI_Bcast(buf, bsize, MPI_CHAR, root, comm);

  // non-root procs write buffer to object
  if (myID != root) {
    // create and open an archive for input
    std::string s(buf, bsize);
    std::stringstream ss(s);
    {
      boost::archive::text_iarchive ia(ss);
      // read class state from archive
      ia >> *t;
    }
  }
}

} /* namespace combigrid */
