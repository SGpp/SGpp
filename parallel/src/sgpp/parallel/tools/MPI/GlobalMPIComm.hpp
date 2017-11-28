// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GLOBALMPICOMM_HPP
#define GLOBALMPICOMM_HPP

// MPI Support
#ifdef USE_MPI
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {
// Hack, remove if possible
extern MPICommunicator* myGlobalMPIComm;
}
}
#endif

#endif /* GLOBALMPICOMM_HPP */
