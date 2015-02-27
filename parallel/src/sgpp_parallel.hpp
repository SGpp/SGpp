// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PARALLEL_HPP
#define PARALLEL_HPP

#ifndef _OPENMP
#error "SGpp parallel module requires OpenMP support"
#endif

#include "sgpp/parallel/datadriven/application/LearnerVectorizedIdentity.hpp"
#include "sgpp/parallel/datadriven/application/LearnerVectorizedIdentitySP.hpp"

#include "sgpp/parallel/operation/ParallelOpFactory.hpp"
#include "sgpp/parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp"

#endif /* PARALLEL_HPP */