// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPITaskScheduler.hpp>

namespace sgpp {
namespace datadriven {
void MPITaskScheduler::setLearnerInstance(LearnerSGDEOnOffParallel *instance) {
  learnerInstance = instance;
}
}
}
