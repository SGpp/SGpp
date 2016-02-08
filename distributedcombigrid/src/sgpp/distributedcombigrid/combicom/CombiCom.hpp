/*
 * CombiCom.hpp
 *
 *  Created on: Aug 16, 2014
 *      Author: heenemo
 */

#ifndef COMBICOM_HPP_
#define COMBICOM_HPP_

#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/SGrid.hpp"
#include "sgpp/distributedcombigrid/utils/StatsContainer.hpp"
#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGridNonUniform.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGrid.hpp"

namespace {
struct mycomplex {
  combigrid::real r;
  combigrid::real i;
};
}

namespace combigrid {

/*
 template <typename FG_ELEMENT>
 class SGrid;

 template <typename FG_ELEMENT>
 class FullGrid;
 */

class CombiCom {
public:
  // after SGReduce sg will have full size
  // sg will be available on ALL members of comm
  template<typename FG_ELEMENT>
  static void SGReduce(SGrid<FG_ELEMENT>& sg, MPI_Comm comm);

  // reduced fg will ONLY be available on r
  template<typename FG_ELEMENT>
  static void FGReduce(FullGrid<FG_ELEMENT>& fg, RankType r, MPI_Comm comm);

  // reduced fg will be available on all member of comm
  template<typename FG_ELEMENT>
  static void FGAllreduce(FullGrid<FG_ELEMENT>& fg, MPI_Comm comm);

  // multiply dfg with coeff and add to dsg. dfg will not be changed
  template<typename FG_ELEMENT>
  static void distributedLocalReduce(DistributedFullGrid<FG_ELEMENT>& dfg,
      DistributedSparseGrid<FG_ELEMENT>& dsg, real coeff);

  template<typename FG_ELEMENT>
  static void distributedLocalReduceNB(DistributedFullGrid<FG_ELEMENT>& dfg,
      DistributedSparseGrid<FG_ELEMENT>& dsg, real coeff);

  template<typename FG_ELEMENT>
  static void distributedLocalReduceBlock(
      DistributedFullGridNonUniform<FG_ELEMENT>& dfg,
      DistributedSparseGrid<FG_ELEMENT>& dsg, real coeff);

  template<typename FG_ELEMENT>
  static void distributedLocalReduceRed(DistributedFullGrid<FG_ELEMENT>& dfg,
      DistributedSparseGrid<FG_ELEMENT>& dsg, real coeff);

  template<typename FG_ELEMENT>
  static void distributedLocalReduceSGR(DistributedFullGrid<FG_ELEMENT>& dfg,
      DistributedSparseGrid<FG_ELEMENT>& dsg, real coeff);

  // extract subspaces of dfg from
  template<typename FG_ELEMENT>
  static void distributedLocalScatter(DistributedFullGrid<FG_ELEMENT>& dfg,
      DistributedSparseGrid<FG_ELEMENT>& dsg);

  template<typename FG_ELEMENT>
  static void distributedGlobalReduce(DistributedSparseGrid<FG_ELEMENT>& dsg);

  template<typename FG_ELEMENT>
  static void
  distributedGlobalReduce(DistributedSparseGridUniform<FG_ELEMENT>& dsg);

  /* create the global communicators for the global reduce.
   * all processes which have local rank i in their process group will be grouped in a
   * distinct communicator.
   * this will only work if all pgroups have same size and dsg uses the same assignment
   * of subspaces to processes in each pgroup
   */
  // todo: this is just a preliminary function. find a more robust way, by only using
  // information which is stored in a unique global object like theMPISystem
  static void registerDistributedGlobalReduceCommmunicator(int pgroupSize);

private:
  static MPI_Comm getDistributedGlobalReduceCommunicator();
};

template<>
inline void CombiCom::SGReduce<double>(SGrid<double>& sg, MPI_Comm comm) {
  // init all empty subspaces
  for (size_t i = 0; i < sg.getSize(); ++i)
    if (sg.getDataSize(i) == 0)
      sg.initHierarchicalSpace(i, 0.0);

  // erzeuge buffer in der größe vom vollen sg
  std::vector<double> buf(sg.getCombinedDataSize());

  // kopiere werte in buffer an richtige stelle
  size_t idx(0);
  for (size_t i = 0; i < sg.getSize(); ++i) {
    real* data = sg.getData(i);
    for (size_t j = 0; j < sg.getDataSize(i); ++j) {
      buf[idx] = data[j];
      ++idx;
    }
  }

  // mpi allreduce
  MPI_Allreduce( MPI_IN_PLACE, &buf[0], static_cast<int>(buf.size()),
      MPI_DOUBLE,
      MPI_SUM, comm);

  // kopiere buffer in richtige stelle an sg_tmp
  idx = 0;
  for (size_t i = 0; i < sg.getSize(); ++i) {
    for (size_t j = 0; j < sg.getDataSize(i); ++j) {
      sg.getData(i)[j] = buf[idx];
      ++idx;
    }
  }
}

template<>
inline void CombiCom::SGReduce<float>(SGrid<float>& sg, MPI_Comm comm) {

  // init all empty subspaces
  for (size_t i = 0; i < sg.getSize(); ++i)
    if (sg.getDataSize(i) == 0)
      sg.initHierarchicalSpace(i, 0.0);

  // erzeuge buffer in der größe vom vollen sg
  std::vector<float> buf(sg.getCombinedDataSize());

  // kopiere werte in buffer an richtige stelle
  size_t idx(0);
  for (size_t i = 0; i < sg.getSize(); ++i) {
    float* data = sg.getData(i);
    for (size_t j = 0; j < sg.getDataSize(i); ++j) {
      buf[idx] = data[j];
      ++idx;
    }
  }

  // mpi allreduce
  MPI_Allreduce( MPI_IN_PLACE, &buf[0], static_cast<int>(buf.size()), MPI_FLOAT,
  MPI_SUM, comm);

  // kopiere buffer in richtige stelle an sg_tmp
  idx = 0;
  for (size_t i = 0; i < sg.getSize(); ++i) {
    for (size_t j = 0; j < sg.getDataSize(i); ++j) {
      sg.getData(i)[j] = buf[idx];
      ++idx;
    }
  }
}

template<>
inline void CombiCom::SGReduce<std::complex<double> >(
    SGrid<std::complex<double> >& sg, MPI_Comm comm) {

  int rank;
  MPI_Comm_rank(comm, &rank);
  std::cout << "rank " << rank << " starting SGreduce" << std::endl;

  // init all empty subspaces
  for (size_t i = 0; i < sg.getSize(); ++i)
    if (sg.getDataSize(i) == 0)
      sg.initHierarchicalSpace(i, 0.0);

  // erzeuge buffer in der größe vom vollen sg
  std::vector<mycomplex> buf(sg.getCombinedDataSize());

  // kopiere werte in buffer an richtige stelle
  size_t idx(0);
  for (size_t i = 0; i < sg.getSize(); ++i) {
    for (size_t j = 0; j < sg.getDataSize(i); ++j) {
      buf[idx].r = sg.getData(i)[j].real();
      buf[idx].i = sg.getData(i)[j].imag();
      ++idx;
    }
  }

  // mpi allreduce
  MPI_Allreduce( MPI_IN_PLACE, &buf[0], static_cast<int>(buf.size()),
  MPI_DOUBLE_COMPLEX, MPI_SUM, comm);

  // kopiere buffer in richtige stelle an sg_tmp
  idx = 0;
  for (size_t i = 0; i < sg.getSize(); ++i) {
    for (size_t j = 0; j < sg.getDataSize(i); ++j) {
      sg.getData(i)[j].real(buf[idx].r);
      sg.getData(i)[j].imag(buf[idx].i);
      ++idx;
    }
  }
}

template<>
inline void CombiCom::FGReduce<double>(FullGrid<double>& fg, RankType r,
    MPI_Comm comm) {
  if (!fg.isGridCreated())
    fg.createFullGrid();

  std::vector<double>& sendbuf = fg.getElementVector();
  std::vector<double> recvbuf(sendbuf.size());

  int myrank;
  MPI_Comm_rank(comm, &myrank);

  if (myrank == r) {
    MPI_Reduce( MPI_IN_PLACE, &sendbuf[0], static_cast<int>(sendbuf.size()),
    MPI_DOUBLE, MPI_SUM, r, comm);
  } else {
    MPI_Reduce(&sendbuf[0], &recvbuf[0], static_cast<int>(sendbuf.size()),
    MPI_DOUBLE, MPI_SUM, r, comm);
  }
}

template<>
inline void CombiCom::FGAllreduce<double>(FullGrid<double>& fg, MPI_Comm comm) {
  if (!fg.isGridCreated())
    fg.createFullGrid();

  std::vector<double>& sendbuf = fg.getElementVector();
  std::vector<double> recvbuf(sendbuf.size());

  int myrank;
  MPI_Comm_rank(comm, &myrank);

  MPI_Allreduce( MPI_IN_PLACE, &sendbuf[0], static_cast<int>(sendbuf.size()),
  MPI_DOUBLE, MPI_SUM, comm);
}

template<>
inline void CombiCom::FGReduce<std::complex<double> >(
    FullGrid<std::complex<double> >& fg, RankType r, MPI_Comm comm) {
  if (!fg.isGridCreated())
    fg.createFullGrid();

  std::vector<complex>& sendbuf = fg.getElementVector();

  std::vector<complex> recvbuf(sendbuf.size());

  int myrank;
  MPI_Comm_rank(comm, &myrank);

  if (myrank == r) {
    MPI_Reduce( MPI_IN_PLACE, &sendbuf[0], static_cast<int>(sendbuf.size()),
    MPI_DOUBLE_COMPLEX, MPI_SUM, r, comm);
  } else {
    MPI_Reduce(&sendbuf[0], &recvbuf[0], static_cast<int>(sendbuf.size()),
    MPI_DOUBLE_COMPLEX, MPI_SUM, r, comm);
  }
}

template<>
inline void CombiCom::FGAllreduce<std::complex<double> >(
    FullGrid<std::complex<double> >& fg, MPI_Comm comm) {
  if (!fg.isGridCreated())
    fg.createFullGrid();

  std::vector<complex>& sendbuf = fg.getElementVector();
  std::vector<complex> recvbuf(sendbuf.size());

  int myrank;
  MPI_Comm_rank(comm, &myrank);

  MPI_Allreduce( MPI_IN_PLACE, &sendbuf[0], static_cast<int>(sendbuf.size()),
  MPI_DOUBLE_COMPLEX, MPI_SUM, comm);
}

template<typename FG_ELEMENT>
void CombiCom::SGReduce(SGrid<FG_ELEMENT>& sg, MPI_Comm comm) {
  assert(!"this type is not yet implemented");
}

template<typename FG_ELEMENT>
void CombiCom::FGReduce(FullGrid<FG_ELEMENT>& fg, RankType r, MPI_Comm comm) {
  assert(!"this type is not yet implemented");
}

template<typename FG_ELEMENT>
void CombiCom::FGAllreduce(FullGrid<FG_ELEMENT>& fg, MPI_Comm comm) {
  assert(!"this type is not yet implemented");
}

template<typename FG_ELEMENT>
void CombiCom::distributedLocalReduce(DistributedFullGrid<FG_ELEMENT>& dfg,
    DistributedSparseGrid<FG_ELEMENT>& dsg, real coeff) {
  assert(dfg.getDimension() == dsg.getDim());
  assert(dfg.returnBoundaryFlags() == dsg.getBoundaryVector());
  assert(dfg.getCommunicatorSize() == dsg.getCommunicatorSize());
  /* todo: check if communicators are equal in the sense that they are a set of the
   * same processes. when dfg creates the cartcomm from lcomm, it changes the
   * communicator (which is just an int value), so simply comparing the communicators
   * is not possible. an easy workaroun could be to additionally store the original comm
   * which has been used in the constructor
   */

  // copy data into subspace datastructures
  dfg.fillSubspaces();

  // get a list of all the subspaces in dfg
  std::vector<LevelVector> dfgSubspacesLevels;
  dfg.getSubspacesLevelVectors(dfgSubspacesLevels);

  // rank of each process
  int rank = dfg.getMpiRank();

  // loop over all subspaces include in dfg
  // and gather those which are contained in the sparse grid
  for (size_t i = 0; i < dfgSubspacesLevels.size(); ++i) {
    const LevelVector& gatherL = dfgSubspacesLevels[i];

    if (!dsg.isContained(gatherL))
      continue;

    std::vector<FG_ELEMENT> buf;

    // index of the current subspace in dsg
    size_t dsgSubID = dsg.getIndex(gatherL);

    // destination process to store the subspace
    RankType dst = dsg.getRank(dsgSubID);

    // start timing
    double tstart = MPI_Wtime();

    dfg.gatherSubspace(gatherL, dst, buf);

    // stop timing
    double time = MPI_Wtime() - tstart;

    // output min, max, avg
    double min, max, sum;
    MPI_Reduce(&time, &min, 1, MPI_DOUBLE, MPI_MIN, dst, dfg.getCommunicator());
    MPI_Reduce(&time, &max, 1, MPI_DOUBLE, MPI_MAX, dst, dfg.getCommunicator());
    MPI_Reduce(&time, &sum, 1, MPI_DOUBLE, MPI_SUM, dst, dfg.getCommunicator());

    double avg = sum / static_cast<double>(dfg.getCommunicatorSize());

    if (rank == dst) {
      std::cout << gatherL << " size = " << buf.size() << " times = " << min
          << " " << max << " " << avg << std::endl;
    }

    // add the subspace data of dfg to dsg weighted by coeff
    if (rank == dst) {
      // this will only init the subspace if not yet initialized
      dsg.initSubspace(dsgSubID, 0.0);
      FG_ELEMENT* data = dsg.getData(dsgSubID);

      for (size_t j = 0; j < dsg.getDataSize(dsgSubID); ++j) {
        data[j] = coeff * buf[j];
      }
    }
  }
}

template<typename FG_ELEMENT>
void CombiCom::distributedLocalReduceRed(DistributedFullGrid<FG_ELEMENT>& dfg,
    DistributedSparseGrid<FG_ELEMENT>& dsg, real coeff) {
  assert(dfg.getDimension() == dsg.getDim());
  assert(dfg.returnBoundaryFlags() == dsg.getBoundaryVector());
  assert(dfg.getCommunicatorSize() == dsg.getCommunicatorSize());
  /* todo: check if communicators are equal in the sense that they are a set of the
   * same processes. when dfg creates the cartcomm from lcomm, it changes the
   * communicator (which is just an int value), so simply comparing the communicators
   * is not possible. an easy workaroun could be to additionally store the original comm
   * which has been used in the constructor
   */

  // copy data into subspace datastructures
  dfg.fillSubspaces();

  // get a list of all the subspaces in dfg
  std::vector<LevelVector> dfgSubspacesLevels;
  dfg.getSubspacesLevelVectors(dfgSubspacesLevels);

  // rank of each process
  int rank = dfg.getMpiRank();

  // create buf that can store largest subspace
  std::vector<FG_ELEMENT> buf(dfg.getMaxSubspaceSize());

  // loop over all subspaces include in dfg
  // and gather those which are contained in the sparse grid
  for (size_t i = 0; i < dfgSubspacesLevels.size(); ++i) {
    const LevelVector& gatherL = dfgSubspacesLevels[i];

    if (!dsg.isContained(gatherL))
      continue;

    // index of the current subspace in dsg
    size_t dsgSubID = dsg.getIndex(gatherL);

    // destination process to store the subspace
    RankType dst = dsg.getRank(dsgSubID);

    // start timing
    //double tstart = MPI_Wtime();

    dfg.gatherSubspaceRed(gatherL, dst, buf);

    // stop timing
    //double time = MPI_Wtime() - tstart;

    // output min, max, avg
    /*
     double min, max, sum;
     MPI_Reduce( &time, &min, 1, MPI_DOUBLE, MPI_MIN, dst, dfg.getCommunicator() );
     MPI_Reduce( &time, &max, 1, MPI_DOUBLE, MPI_MAX, dst, dfg.getCommunicator() );
     MPI_Reduce( &time, &sum, 1, MPI_DOUBLE, MPI_SUM, dst, dfg.getCommunicator() );
     */

    //double avg = sum / static_cast<double>( dfg.getCommunicatorSize() );
    /*
     if( rank == dst ){
     std::cout << gatherL << " size = " << buf.size()
     << " times = " << min << " " << max << " " << avg << std::endl;
     }*/

    // add the subspace data of dfg to dsg weighted by coeff
    if (rank == dst) {
      // this will only init the subspace if not yet initialized
      dsg.initSubspace(dsgSubID, 0.0);
      FG_ELEMENT* data = dsg.getData(dsgSubID);

      for (size_t j = 0; j < dsg.getDataSize(dsgSubID); ++j) {
        data[j] = coeff * buf[j];
      }
    }
  }
}

template<typename FG_ELEMENT>
void CombiCom::distributedLocalReduceBlock(
    DistributedFullGridNonUniform<FG_ELEMENT>& dfg,
    DistributedSparseGrid<FG_ELEMENT>& dsg, real coeff) {
  theStatsContainer()->setTimerStart("local_reduce_fill_subspaces");
  // copy data into subspace datastructures
  dfg.fillSubspaces();
  theStatsContainer()->setTimerStop("local_reduce_fill_subspaces");

  theStatsContainer()->setTimerStart("local_reduce_calc_dependencies");
  // get a list of all the subspaces in dfg
  std::vector<LevelVector> dfgSubspacesLevels;
  dfg.getSubspacesLevelVectors(dfgSubspacesLevels);

  // buffers for send and recvdata, one for each process
  std::vector<std::vector<FG_ELEMENT> > senddata(dfg.getCommunicatorSize());
  std::vector<std::vector<int> > sendsizes(dfg.getCommunicatorSize());
  std::vector<std::vector<int> > sendsubspaces(dfg.getCommunicatorSize());
  std::vector<std::vector<FG_ELEMENT> > recvdata(dfg.getCommunicatorSize());
  std::vector<std::vector<int> > recvsizes(dfg.getCommunicatorSize());
  std::vector<std::vector<int> > recvsubspaces(dfg.getCommunicatorSize());

  // loop over all subspaces include in dfg
  // and gather those which are contained in the sparse grid
  for (size_t i = 0; i < dfgSubspacesLevels.size(); ++i) {
    const LevelVector& gatherL = dfgSubspacesLevels[i];

    if (!dsg.isContained(gatherL))
      continue;

    std::vector<FG_ELEMENT> buf;

    // index of the current subspace in dsg
    size_t dsgSubID = dsg.getIndex(gatherL);

    // destination process to store the subspace
    RankType dst = dsg.getRank(dsgSubID);

    // get send and recvdata for this subspace
    dfg.gatherSubspaceBlock(gatherL, dst, senddata, sendsizes, sendsubspaces,
        recvsizes, recvsubspaces);
  }
  theStatsContainer()->setTimerStop("local_reduce_calc_dependencies");

  theStatsContainer()->setTimerStart("local_reduce_exchange_data");
  // start send and recv operations to all other processes
  std::vector<MPI_Request> requests;
  size_t totalsendsize(0), totalrecvsize(0);
  for (int r = 0; r < dfg.getCommunicatorSize(); ++r) {
    MPI_Request sendrequest;
    MPI_Isend(senddata[r].data(), int(senddata[r].size()),
    MPI_INT, r, 0, dfg.getCommunicator(), &sendrequest);
    requests.push_back(sendrequest);

    totalsendsize += senddata[r].size();

    // calc recv data size
    int rsize = 0;
    for (auto s : recvsizes[r])
      rsize += s;
    recvdata[r].resize(rsize);

    MPI_Request recvrequest;
    MPI_Irecv(recvdata[r].data(), rsize, MPI_INT, r, 0, dfg.getCommunicator(),
        &recvrequest);
    requests.push_back(recvrequest);

    totalrecvsize += rsize;
  }

  MPI_Waitall(int(requests.size()), &requests[0], MPI_STATUSES_IGNORE);
  theStatsContainer()->setTimerStop("local_reduce_exchange_data");

  int rank = dfg.getMpiRank();
  for (int r = 0; r < dfg.getCommunicatorSize(); ++r) {
    if (r == rank) {
      std::cout << "rank " << rank << " tot send " << totalsendsize
          << " tot recv " << totalrecvsize << std::endl;
    }

    MPI_Barrier(dfg.getCommunicator());
  }

  // todo: for own dsg subspaces -> put togheter subspace data from recvdata
}

template<typename FG_ELEMENT>
void CombiCom::distributedLocalReduceNB(DistributedFullGrid<FG_ELEMENT>& dfg,
    DistributedSparseGrid<FG_ELEMENT>& dsg, real coeff) {
  assert(dfg.getDimension() == dsg.getDim());
  assert(dfg.returnBoundaryFlags() == dsg.getBoundaryVector());
  assert(dfg.getCommunicatorSize() == dsg.getCommunicatorSize());
  /* todo: check if communicators are equal in the sense that they are a set of the
   * same processes. when dfg creates the cartcomm from lcomm, it changes the
   * communicator (which is just an int value), so simply comparing the communicators
   * is not possible. an easy workaroun could be to additionally store the original comm
   * which has been used in the constructor
   */

  theStatsContainer()->setTimerStart("local_reduce_fill_subspaces");
  // copy data into subspace datastructures
  dfg.fillSubspaces();
  theStatsContainer()->setTimerStop("local_reduce_fill_subspaces");

  theStatsContainer()->setTimerStart("local_reduce_exchange_data");

  // get a list of all the subspaces in dfg
  std::vector<LevelVector> dfgSubspacesLevels;
  dfg.getSubspacesLevelVectors(dfgSubspacesLevels);

  // create buffers for all subspaces, but don't init size
  std::vector<std::vector<FG_ELEMENT> > buffers(dfgSubspacesLevels.size());

  // create storage for all the requests
  std::vector<MPI_Request> requests;

  // loop over all subspaces include in dfg
  // and gather those which are contained in the sparse grid
  for (size_t i = 0; i < dfgSubspacesLevels.size(); ++i) {
    const LevelVector& gatherL = dfgSubspacesLevels[i];

    if (!dsg.isContained(gatherL))
      continue;

    // get buffer for this subspace
    std::vector<FG_ELEMENT>& buf = buffers[i];

    // index of the current subspace in dsg
    size_t dsgSubID = dsg.getIndex(gatherL);

    // destination process to store the subspace
    RankType dst = dsg.getRank(dsgSubID);

    // in request only request for this process are returned
    dfg.gatherSubspaceNB(gatherL, dst, buf, requests);
  }

  // wait for requests
  if (requests.size() > 0)
    MPI_Waitall(static_cast<int>(requests.size()), &requests[0],
        MPI_STATUSES_IGNORE);

  theStatsContainer()->setTimerStop("local_reduce_exchange_data");

  theStatsContainer()->setTimerStart("local_reduce_add_buffers");
  // each process adds the buffers it is responsible for to the corresponding
  // subpsace in dsg
  for (size_t i = 0; i < dfgSubspacesLevels.size(); ++i) {
    if (buffers[i].size() == 0)
      continue;

    const LevelVector& gatherL = dfgSubspacesLevels[i];

    // index of the current subspace in dsg
    size_t dsgSubID = dsg.getIndex(gatherL);

    // this will only init the subspace if not yet initialized
    dsg.initSubspace(dsgSubID, 0.0);

    FG_ELEMENT* data = dsg.getData(dsgSubID);

    // get buffer for this subspace
    std::vector<FG_ELEMENT>& buf = buffers[i];

    for (size_t j = 0; j < dsg.getDataSize(dsgSubID); ++j) {
      data[j] = coeff * buf[j];
    }
  }

  theStatsContainer()->setTimerStop("local_reduce_add_buffers");
}

template<typename FG_ELEMENT>
void CombiCom::distributedLocalReduceSGR(DistributedFullGrid<FG_ELEMENT>& dfg,
    DistributedSparseGrid<FG_ELEMENT>& dsg, real coeff) {
  assert(dfg.getDimension() == dsg.getDim());
  assert(dfg.returnBoundaryFlags() == dsg.getBoundaryVector());
  assert(dfg.getCommunicatorSize() == dsg.getCommunicatorSize());
  /* todo: check if communicators are equal in the sense that they are a set of the
   * same processes. when dfg creates the cartcomm from lcomm, it changes the
   * communicator (which is just an int value), so simply comparing the communicators
   * is not possible. an easy workaroun could be to additionally store the original comm
   * which has been used in the constructor
   */

  // get list of subspaces in dfg and sort according to rank
  std::vector<LevelVector> dfgSubspacesLevels;
  dfg.getSubspacesLevelVectors(dfgSubspacesLevels);

  std::vector<LevelVector> commonSubspaces;
  std::vector<int> commonSubspacesRanks;
  std::vector<size_t> commonSubspacesSizes;

  for (size_t i = 0; i < dfgSubspacesLevels.size(); ++i) {
    const LevelVector& gatherL = dfgSubspacesLevels[i];

    if (!dsg.isContained(gatherL))
      continue;

    commonSubspaces.push_back(gatherL);
    commonSubspacesRanks.push_back(dsg.getRank(gatherL));
    commonSubspacesSizes.push_back(dsg.getSubspaceSize(gatherL));
  }

  // create buffer for each rank
  std::vector<std::vector<FG_ELEMENT> > buffers(dfg.getCommunicatorSize());

  // for each rank
  for (int r = 0; r < dfg.getCommunicatorSize(); ++r) {
    std::vector<LevelVector> mySubspaces;
    std::vector<size_t> mySizes;

    size_t bsize = 0;
    for (size_t i = 0; i < commonSubspaces.size(); ++i) {
      if (commonSubspacesRanks[i] == r) {
        mySubspaces.push_back(commonSubspaces[i]);
        mySizes.push_back(commonSubspacesSizes[i]);
        bsize += commonSubspacesSizes[i];
      }
    }

    // resize buffer
    buffers[r].resize(bsize, FG_ELEMENT(0));

    // for each subspace for this rank
    // todo: put local data into right position of buffer

    // mpi reduce on buffer
    if (dfg.getMpiRank() == r) {
      MPI_Reduce( MPI_IN_PLACE, buffers[r].data(), bsize, dfg.getMPIDatatype(),
      MPI_SUM, r, dfg.getCommunicator());
    } else {
      MPI_Reduce(buffers[r].data(), buffers[r].data(), bsize,
          dfg.getMPIDatatype(), MPI_SUM, r, dfg.getCommunicator());
    }

    if (dfg.getMpiRank() != r)
      buffers[r].resize(0);
  }
}

template<typename FG_ELEMENT>
void CombiCom::distributedLocalScatter(DistributedFullGrid<FG_ELEMENT>& dfg,
    DistributedSparseGrid<FG_ELEMENT>& dsg) {
  assert(dfg.getDimension() == dsg.getDim());
  assert(dfg.returnBoundaryFlags() == dsg.getBoundaryVector());
  assert(dfg.getCommunicatorSize() == dsg.getCommunicatorSize());
  /* todo: check if communicators are equal in the sense that they are a set of the
   * same processes. when dfg creates the cartcomm from lcomm, it changes the
   * communicator (which is just an int value), so simply comparing the communicators
   * is not possible. an easy workaroun could be to additionally store the original comm
   * which has been used in the constructor
   */

  theStatsContainer()->setTimerStart("local_scatter_exchange_data");
  // get a list of all the subspaces in dfg
  std::vector<LevelVector> dfgSubspacesLevels;
  dfg.getSubspacesLevelVectors(dfgSubspacesLevels);

  // loop over all subspaces include in dfg
  // and scatter all sparse grid subspaces contained in dfg those which are contained in the sparse grid
  for (size_t i = 0; i < dfgSubspacesLevels.size(); ++i) {
    const LevelVector& scatterL = dfgSubspacesLevels[i];

    if (!dsg.isContained(scatterL))
      continue;

    // index of the current subspace in dsg
    size_t dsgSubID = dsg.getIndex(scatterL);

    // process that stores the subspace
    RankType src = dsg.getRank(dsgSubID);

    const std::vector<real>& buf = dsg.getDataVector(dsgSubID);

    dfg.scatterSubspace(scatterL, src, buf);
  }
  theStatsContainer()->setTimerStop("local_scatter_exchange_data");

  // copy data back into dfg
  theStatsContainer()->setTimerStart("local_scatter_write_back");
  dfg.writeBackSubspaces();
  theStatsContainer()->setTimerStop("local_scatter_write_back");

  // clear subspace containers and release memory
  dfg.clearSubspaces();
}

template<typename FG_ELEMENT>
void CombiCom::distributedGlobalReduce(DistributedSparseGrid<FG_ELEMENT>& dsg) {
  // first time -> register communicators
  CombiCom::registerDistributedGlobalReduceCommmunicator(
      dsg.getCommunicatorSize());

  // get rank in pgroup communicator.
  int lrank;
  MPI_Comm_rank(dsg.getCommunicator(), &lrank);

  std::vector<MPI_Request> myrequests;

  for (size_t i = 0; i < dsg.getNumSubspaces(); ++i) {
    // skip all subspaces which are not stored on lrank
    if (dsg.getRank(i) != lrank)
      continue;

    // get global communicator for this subspace
    MPI_Comm mycomm = CombiCom::getDistributedGlobalReduceCommunicator();

    // make sure that subspace is initialized. not all subspaces will be initialized
    // after local reduce. this will not overwrite an already initialized subspace
    dsg.initSubspace(i, 0.0);

    FG_ELEMENT* buf = dsg.getData(i);
    int bsize = int(dsg.getDataSize(i));
    MPI_Datatype dtype = abstraction::getMPIDatatype(
        abstraction::getabstractionDataType<FG_ELEMENT>());

    // mpi allreduce
    if (USE_NONBLOCKING_MPI_COLLECTIVE) {
      MPI_Request request;
      MPI_Iallreduce( MPI_IN_PLACE, buf, bsize, dtype, MPI_SUM, mycomm,
          &request);
      myrequests.push_back(request);
    } else {
      MPI_Allreduce( MPI_IN_PLACE, buf, bsize, dtype, MPI_SUM, mycomm);
    }
  }

  if (USE_NONBLOCKING_MPI_COLLECTIVE) {
    MPI_Waitall(int(myrequests.size()), &myrequests[0], MPI_STATUS_IGNORE);
  }
}

template<typename FG_ELEMENT>
void CombiCom::distributedGlobalReduce(
    DistributedSparseGridUniform<FG_ELEMENT>& dsg) {
  // get global communicator for this operation
  MPI_Comm mycomm = CombiCom::getDistributedGlobalReduceCommunicator();

  assert(mycomm != MPI_PROC_NULL);

  /* get sizes of all partial subspaces in communicator
   * we have to do this, because size information of uninitialized subspaces
   * is not available in dsg. at the moment this information is only available
   * in dfg.
   */
  std::vector<int> subspaceSizes(dsg.getNumSubspaces());
  for (size_t i = 0; i < subspaceSizes.size(); ++i) {
    // MPI does not have a real size_t equivalent. int should work in most cases
    // if not we can at least detect this with an assert
    assert(dsg.getDataSize(i) <= INT_MAX);

    subspaceSizes[i] = int(dsg.getDataSize(i));
  }

  MPI_Allreduce( MPI_IN_PLACE, subspaceSizes.data(), int(subspaceSizes.size()),
  MPI_INT, MPI_MAX, mycomm);

  // check for implementation errors, the reduced subspace size should not be
  // different from the size of already initialized subspaces
  int bsize = 0;
  for (size_t i = 0; i < subspaceSizes.size(); ++i) {
    bool check = (subspaceSizes[i] == 0 || dsg.getDataSize(i) == 0
        || subspaceSizes[i] == int(dsg.getDataSize(i)));
    if (!check) {
      int rank;
      MPI_Comm_rank( MPI_COMM_WORLD, &rank);
      std::cout << "l = " << dsg.getLevelVector(i) << " " << "rank = " << rank
          << " " << "ssize = " << subspaceSizes[i] << " " << "dsize = "
          << dsg.getDataSize(i) << std::endl;
      assert(false);
    }
    bsize += subspaceSizes[i];
  }

  // put subspace data into buffer
  std::vector<FG_ELEMENT> buf(bsize, FG_ELEMENT(0));
  {
    typename std::vector<FG_ELEMENT>::iterator buf_it = buf.begin();
    for (size_t i = 0; i < dsg.getNumSubspaces(); ++i) {
      std::vector<FG_ELEMENT>& subspaceData = dsg.getDataVector(i);

      // if subspace does not exist on this process this part of the buffer is
      // left empty
      if (subspaceData.size() == 0) {
        buf_it += subspaceSizes[i];
      }
      for (size_t j = 0; j < subspaceData.size(); ++j) {
        *buf_it = subspaceData[j];
        ++buf_it;
      }
    }
  }

  MPI_Datatype dtype = abstraction::getMPIDatatype(
      abstraction::getabstractionDataType<FG_ELEMENT>());
  MPI_Allreduce( MPI_IN_PLACE, buf.data(), bsize, dtype, MPI_SUM, mycomm);

  // extract subspace data
  {
    typename std::vector<FG_ELEMENT>::iterator buf_it = buf.begin();
    for (size_t i = 0; i < dsg.getNumSubspaces(); ++i) {
      std::vector<FG_ELEMENT>& subspaceData = dsg.getDataVector(i);

      // if subspace does not exist on this process this part of the buffer is
      // skipped
      if (subspaceData.size() == 0) {
        buf_it += subspaceSizes[i];
      }
      for (size_t j = 0; j < subspaceData.size(); ++j) {
        subspaceData[j] = *buf_it;
        ++buf_it;
      }
    }
  }
}

} /* namespace combigrid */

#endif /* COMBICOM_HPP_ */
