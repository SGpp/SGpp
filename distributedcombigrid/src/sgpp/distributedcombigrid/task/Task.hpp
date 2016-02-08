/*
 * Task.hpp
 *
 *  Created on: Jul 10, 2014
 *      Author: heenemo
 */

#ifndef TASK_HPP_
#define TASK_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <string>
#include <vector>
#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LoadModel.hpp"

namespace combigrid {

/** 
 * Interface for tasks
 */
class Task {
protected:
  Task();

  Task(DimType dim, LevelVector& l, LevelVector& nmax, LevelVector& lmin,
      std::vector<bool>& boundary, real coeff, LoadModel* loadModel);
public:
  virtual ~Task();
  // pgroup receive task
  static void receive(Task** t, RankType source, CommunicatorType& comm);

  // manager send task to pgroup root
  static void send(Task** t, RankType dest, CommunicatorType& comm);

  //broadcast task
  static void broadcast(Task** t, RankType root, CommunicatorType& comm);

  inline DimType getDim() const;

  inline const LevelVector& getLevelVector() const;

  inline const std::vector<bool>& getBoundary() const;

  inline int getID();

  virtual void run(CommunicatorType& lcomm) = 0;

  virtual void init(CommunicatorType& lcomm) = 0;

  inline real estimateRuntime() const;

  inline bool isFinished() const;

  inline void setFinished(bool finished);

  virtual void getFullGrid(FullGrid<CombiDataType>& fg, RankType lroot,
      CommunicatorType& lcomm) = 0;

  virtual DistributedFullGrid<CombiDataType>& getDistributedFullGrid() = 0;

private:
  friend class boost::serialization::access;

  // serialize
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);

protected:
  DimType dim_;

  LevelVector l_;       // levelvector of partial solution

  std::vector<bool> boundary_;

  int id_;              // unique id of task, same on manager and worker

  static int count;

  LoadModel* loadModel_;

  bool isFinished_;
};

typedef std::vector<Task*> TaskContainer;

template<class Archive>
void Task::serialize(Archive & ar, const unsigned int version) {
  ar & dim_;
  ar & id_;
  ar & l_;
  ar & boundary_;
}

inline DimType Task::getDim() const {
  return dim_;
}

inline const LevelVector& Task::getLevelVector() const {
  return l_;
}

inline const std::vector<bool>& Task::getBoundary() const {
  return boundary_;
}

inline int Task::getID() {
  return id_;
}

inline bool Task::isFinished() const {
  return isFinished_;
}

inline void Task::setFinished(bool finished) {
  isFinished_ = finished;
}

inline real Task::estimateRuntime() const {
  return loadModel_->eval(l_);
}

} /* namespace combigrid */

#endif /* TASK_HPP_ */
