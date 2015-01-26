/* ****************************************************************************
* Copyright (C) 201 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef GLOBALMPICOMM_HPP
#define GLOBALMPICOMM_HPP

// MPI Support
#ifdef USE_MPI
namespace sg {
  namespace parallel {
    // @todo MPI (heinecke) Hack, remove if possible
    extern MPICommunicator* myGlobalMPIComm;
  }
}
#endif


#endif /* GLOBALMPICOMM_HPP */
