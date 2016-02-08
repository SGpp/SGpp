/*
 * StatsContainer.hpp
 *
 *  Created on: Jul 9, 2015
 *      Author: heenemo
 */

#ifndef STATSCONTAINER_HPP_
#define STATSCONTAINER_HPP_

#include <map>
#include <chrono>
#include <cassert>
#include <string>
#include <fstream>
#include <iostream>

#include "sgpp/distributedcombigrid/utils/StatsContainerID.hpp"

#define TIMING

namespace combigrid {

/* a class for timing and storing stats
 *
 * times will be measured in seconds
 *
 * output order will be newest entries first
 *
 * values can be overwritten
 *
 * starting or stopping a timer twice will lead to error
 *
 * also trying to get a duration of a non-existing timer will lead to error
 */
class StatsContainer {

 public:
  inline void setTimerStart(const std::string& timerName);

  inline void setTimerStop(const std::string& timerName);

  inline double getDuration(const std::string& timerName) const;

  inline void setValue(const std::string& valueName, double value);

  inline double getValue(const std::string& valueName) const;

  virtual ~StatsContainer() {

  };

  inline void save(const std::string& filename) const;

 private:
  friend StatsContainerID theStatsContainer();

  explicit StatsContainer() {

  }

  std::map<std::string, std::chrono::high_resolution_clock::time_point>
  startTimes_;

  std::map<std::string, std::chrono::high_resolution_clock::time_point>
  stopTimes_;

  std::map<std::string, double> durations_;

  std::map<std::string, double> stats_;
};

#ifdef TIMING
inline void StatsContainer::setTimerStart(const std::string& timerName) {
  using namespace std::chrono;
  // check if timer has already been started
  std::map<std::string, high_resolution_clock::time_point>::const_iterator found =
    startTimes_.find(timerName);

  //assert( found == startTimes_.end() && "timer already started" );
  if (!(found == startTimes_.end())) {
    std::cout << "Warning: timer " << timerName << " already started"
              << std::endl;
  }

  startTimes_[timerName] = high_resolution_clock::now();
}

inline void StatsContainer::setTimerStop(const std::string& timerName) {
  using namespace std::chrono;
  // check if timer has been started
  std::map<std::string, high_resolution_clock::time_point>::const_iterator found =
    startTimes_.find(timerName);

  //assert( found != startTimes_.end() && "timer has not been started" );
  if (found == startTimes_.end()) {
    std::cout << "Warning: timer " << timerName << "  has not been started"
              << std::endl;
  }

  // check if timer has already been stopped
  std::map<std::string, high_resolution_clock::time_point>::const_iterator found2
    =
      stopTimes_.find(timerName);

  //assert( found2 == stopTimes_.end() && "timer was already stopped" );
  if (found2 != stopTimes_.end()) {
    std::cout << "Warning: timer " << timerName << " timer was already stopped"
              << std::endl;
  }

  high_resolution_clock::time_point stopTime = high_resolution_clock::now();
  stopTimes_[timerName] = stopTime;

  // calc duration
  duration<double> dur = stopTime - startTimes_[timerName];
  durations_[timerName] = dur.count();
}

inline double StatsContainer::getDuration(const std::string& timerName) const {
  using namespace std::chrono;
  // check if timer exists
  // check if timer has already been stopped
  std::map<std::string, double>::const_iterator found = durations_.find(
        timerName);

  if (found == durations_.end()) {
    std::cout << "Error: value " << timerName << " does not exist!"
              << std::endl;
    return -1;
  }

  return durations_.at(timerName);
}

inline void StatsContainer::setValue(const std::string& valueName,
                                     double value) {
  // check if value exists
  std::map<std::string, double>::const_iterator found = stats_.find(valueName);
  //assert( found == stats_.end() && "value already exists" );

  stats_[valueName] = value;
}

inline double StatsContainer::getValue(const std::string& valueName) const {
  // check if value exists
  std::map<std::string, double>::const_iterator found = stats_.find(valueName);

  if (found == stats_.end()) {
    std::cout << "Error: value " << valueName << " does not exist!"
              << std::endl;
    return -1;
  }

  return stats_.at(valueName);
}

inline void StatsContainer::save(const std::string& filename) const {
  std::ofstream ofs(filename.c_str());

  for (std::map<std::string, double>::const_iterator it = durations_.begin();
       it != durations_.end(); ++it) {
    ofs << it->first << "\t" << it->second << std::endl;
  }

  for (std::map<std::string, double>::const_iterator it = stats_.begin();
       it != stats_.end(); ++it) {
    ofs << it->first << "\t" << it->second << std::endl;
  }

  ofs.close();
}

/* returns a handle to the (global) StatsContainer. when always using this
 * handle only one StatsContainer
 */
inline StatsContainerID theStatsContainer() {
  static StatsContainerID stats(new StatsContainer());
  return stats;
}
#else
inline void StatsContainer::setTimerStart( const std::string& timerName ) {}

inline void StatsContainer::setTimerStop( const std::string& timerName ) {}

inline double StatsContainer::getDuration( const std::string& timerName )
const {
  return 0.0;
}

inline void StatsContainer::setValue( const std::string& valueName,
                                      double value ) {}

inline double StatsContainer::getValue( const std::string& valueName) const {
  return 0.0;
}

inline void StatsContainer::save( const std::string& filename ) const {}

/* returns a handle to the (global) StatsContainer. when always using this
 * handle only one StatsContainer
 */
inline StatsContainerID theStatsContainer() {
  static StatsContainerID stats( new StatsContainer() );
  return stats;
}
#endif

}
//end namespace combigrid

#endif /* STATSCONTAINER_HPP_ */
