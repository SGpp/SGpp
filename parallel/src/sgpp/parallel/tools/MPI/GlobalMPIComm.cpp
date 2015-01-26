/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/parallel/tools/MPI/MPICommunicator.hpp>

namespace sg {
  namespace parallel {
    // @todo MPI (heinecke) Hack, remove if possible
    MPICommunicator* myGlobalMPIComm;
  }
}
