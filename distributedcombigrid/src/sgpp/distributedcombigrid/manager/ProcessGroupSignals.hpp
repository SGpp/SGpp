/*
 * ProcessGroupCommands.hpp
 *
 *  Created on: Jun 24, 2014
 *      Author: heenemo
 */

#ifndef PROCESSGROUPSIGNALS_HPP_
#define PROCESSGROUPSIGNALS_HPP_

namespace combigrid {

// attention: changing SignalType might require changing the MPI Type as well
typedef int SignalType;

const SignalType RUN_FIRST = 0;
const SignalType RUN_NEXT = 1;
const SignalType EVAL = 2;
const SignalType GRID_EVAL = 3;
const SignalType COMBINE = 4;
const SignalType EXIT = 5;
const SignalType SYNC_TASKS = 6;
const SignalType TEST_CONVERSION = 7;
const SignalType COMBINE_FG = 8;
const SignalType EV_CALC_FG = 9;
const SignalType EV_CALC_FG_INIT = 10;
const SignalType GRID_GATHER = 11;
const SignalType UPDATE_COMBI_PARAMETERS = 12;
const SignalType ADD_TASK = 13;
const SignalType RECOMPUTE = 14;

typedef int NormalizationType;
const NormalizationType NO_NORMALIZATION = 0;
const NormalizationType L1_NORMALIZATION = 1;
const NormalizationType L2_NORMALIZATION = 2;
const NormalizationType EV_NORMALIZATION = 3;

enum TagType {
  signalTag = 0, statusTag = 1
};

// attention: changing StatusType might require changing the MPI Type
typedef int StatusType;

const StatusType PROCESS_GROUP_WAIT = 0;
const StatusType PROCESS_GROUP_BUSY = 1;

} /* namespace combigrid */

#endif /* PROCESSGROUPCOMMANDS_HPP_ */
