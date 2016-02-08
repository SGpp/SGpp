/*
 * Types.hpp
 *
 *  Created on: May 14, 2013
 *      Author: heenemo
 */
#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <cstddef> //contains size_t
#include <complex>
#include <mpi.h>
#include <stdexcept>

namespace combigrid {

typedef double real;

// IndexType must be signed! because we use -1 in some functions as a return
// value. and large enough so that all grid points can be
// numbered. i.e levelsum (with boundary) < 30 for int32 and < 62 for int64
typedef int64_t IndexType;

typedef IndexType LevelType;

typedef std::complex<real> complex;

typedef size_t DimType;

typedef MPI_Comm CommunicatorType;

typedef int RankType;

/* set the datatype for the values stored in any type of grid used in this
 * framework
 */
//typedef complex CombiDataType;
typedef real CombiDataType;

// global defines. we go the more C++ way here and use const variables

/* nonblocking mpi collective calls (MPI_Iallreduce and the likes) usually yield
 * better performance. if you observe problems with these functions uncomment to
 * fall back to blocking counterpart of the function
 */
const bool USE_NONBLOCKING_MPI_COLLECTIVE(true);
}

namespace abstraction {

typedef enum {
  type_unknown = 0,
  type_float,
  type_double,
  type_double_complex,
  type_float_complex
} DataType;

template<class T>
DataType getabstractionDataType() {
  throw new std::invalid_argument("Datatype is not supported!");
}

template<>
inline DataType getabstractionDataType<float>() {
  return abstraction::type_float;
}

template<>
inline DataType getabstractionDataType<double>() {
  return abstraction::type_double;
}

template<>
inline DataType getabstractionDataType<std::complex<double> >() {
  return abstraction::type_double_complex;
}

template<>
inline DataType getabstractionDataType<std::complex<float> >() {
  return abstraction::type_float_complex;
}

static MPI_Datatype getMPIDatatype(abstraction::DataType type) {
  switch (type) {
    case abstraction::type_float:
      return MPI_FLOAT;

    case abstraction::type_double:
      return MPI_DOUBLE;

    case abstraction::type_double_complex:
      return MPI_DOUBLE_COMPLEX;

    case abstraction::type_float_complex:
      return MPI_COMPLEX;

    case abstraction::type_unknown:
      throw new std::invalid_argument("Type unknown ConvertType!");
  };

  throw new std::invalid_argument(
    "MPI_Datatype Convert(abstraction::DataType) failed!");
}

}

#endif /* TYPES_HPP_ */
