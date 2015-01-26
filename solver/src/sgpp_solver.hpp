/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de

#ifndef SOLVER_HPP
#define SOLVER_HPP


#include "sgpp/solver/sle/ConjugateGradients.hpp"
#include "sgpp/solver/sle/BiCGStab.hpp"
#include "sgpp/solver/ode/Euler.hpp"
#include "sgpp/solver/ode/CrankNicolson.hpp"
#include "sgpp/solver/ode/AdamsBashforth.hpp"
#include "sgpp/solver/ode/VarTimestep.hpp"
#include "sgpp/solver/ode/StepsizeControl.hpp"
#include "sgpp/solver/ode/StepsizeControlEJ.hpp"
#include "sgpp/solver/ode/StepsizeControlH.hpp"
#include "sgpp/solver/ode/StepsizeControlMC.hpp"
#include "sgpp/solver/ode/StepsizeControlBDF.hpp"
#include "sgpp/solver/TypesSolver.hpp"

#endif /* SOLVER_HPP */
