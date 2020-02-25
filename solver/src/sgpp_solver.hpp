// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/solver/ode/Euler.hpp>
#include <sgpp/solver/ode/CrankNicolson.hpp>
#include <sgpp/solver/ode/AdamsBashforth.hpp>
#include <sgpp/solver/ode/VarTimestep.hpp>
#include <sgpp/solver/ode/StepsizeControl.hpp>
#include <sgpp/solver/ode/StepsizeControlEJ.hpp>
#include <sgpp/solver/ode/StepsizeControlH.hpp>
#include <sgpp/solver/ode/StepsizeControlMC.hpp>
#include <sgpp/solver/ode/StepsizeControlBDF.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/solver/SLESolverTypeParser.hpp>
#include <sgpp/solver/operation/hash/OperationParabolicPDESolverSystem.hpp>

#endif /* SOLVER_HPP */
