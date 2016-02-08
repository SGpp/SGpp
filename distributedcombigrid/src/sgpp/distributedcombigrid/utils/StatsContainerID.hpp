/*
 * StatsContainerID.hpp
 *
 *  Created on: Jul 9, 2015
 *      Author: heenemo
 */

#ifndef STATSCONTAINERID_HPP_
#define STATSCONTAINERID_HPP_

#include <memory>

namespace combigrid {

class StatsContainer;

typedef std::shared_ptr<StatsContainer> StatsContainerID;

} // namespace combigrid

#endif /* STATSCONTAINERID_HPP_ */
