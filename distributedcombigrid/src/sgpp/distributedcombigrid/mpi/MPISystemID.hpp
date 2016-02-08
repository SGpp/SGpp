/*
 * MPISystemID.hpp
 *
 *  Created on: Jan 23, 2013
 *      Author: mh
 *
 *  Partially copied from the pe Physics Engine class MPISystemID
 */

#ifndef MPISYSTEMID_HPP
#define MPISYSTEMID_HPP

#include <tr1/memory>

namespace combigrid {

class MPISystem;

/*!\brief Handle for the MPI communication system.
 // \ingroup mpi
 */
typedef std::tr1::shared_ptr<MPISystem> MPISystemID;

/*!\brief Handle for the constant MPI communication system.
 // \ingroup mpi
 */
typedef std::tr1::shared_ptr<const MPISystem> ConstMPISystemID;

} // namespace combigrid

#endif // MPISYSTEMID_HPP
