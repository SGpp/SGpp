// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "FitterConfiguration.hpp"

#include <string>

namespace sgpp {
namespace datadriven {

const base::RegularGridConfiguration& FitterConfiguration::getGridConfig() const {
  return gridConfig;
}

const base::AdpativityConfiguration& FitterConfiguration::getRefinementConfig() const {
  return adaptivityConfig;
}

const datadriven::DecompositionConfiguration& FitterConfiguration::getDecompositionConfig() const {
  return decompositionConfig;
}

const solver::SLESolverConfiguration& FitterConfiguration::getSolverRefineConfig() const {
  return solverRefineConfig;
}

const solver::SLESolverConfiguration& FitterConfiguration::getSolverFinalConfig() const {
  return solverFinalConfig;
}

const datadriven::RegularizationConfiguration& FitterConfiguration::getRegularizationConfig()
    const {
  return regularizationConfig;
}

const datadriven::OperationMultipleEvalConfiguration& FitterConfiguration::getMultipleEvalConfig()
    const {
  return multipleEvalConfig;
}

base::RegularGridConfiguration& FitterConfiguration::getGridConfig() {
  return const_cast<base::RegularGridConfiguration&>(
      static_cast<const FitterConfiguration&>(*this).getGridConfig());
}

base::AdpativityConfiguration& FitterConfiguration::getRefinementConfig() {
  return const_cast<base::AdpativityConfiguration&>(
      static_cast<const FitterConfiguration&>(*this).getRefinementConfig());
}

datadriven::DecompositionConfiguration& FitterConfiguration::getDecompositionConfig() {
  return const_cast<datadriven::DecompositionConfiguration&>(
      static_cast<const FitterConfiguration&>(*this).getDecompositionConfig());
}

solver::SLESolverConfiguration& FitterConfiguration::getSolverRefineConfig() {
  return const_cast<solver::SLESolverConfiguration&>(
      static_cast<const FitterConfiguration&>(*this).getSolverRefineConfig());
}

solver::SLESolverConfiguration& FitterConfiguration::getSolverFinalConfig() {
  return const_cast<solver::SLESolverConfiguration&>(
      static_cast<const FitterConfiguration&>(*this).getSolverFinalConfig());
}

datadriven::RegularizationConfiguration& FitterConfiguration::getRegularizationConfig() {
  return const_cast<datadriven::RegularizationConfiguration&>(
      static_cast<const FitterConfiguration&>(*this).getRegularizationConfig());
}

datadriven::OperationMultipleEvalConfiguration& FitterConfiguration::getMultipleEvalConfig() {
  return const_cast<datadriven::OperationMultipleEvalConfiguration&>(
      static_cast<const FitterConfiguration&>(*this).getMultipleEvalConfig());
}
}  // namespace datadriven
}  // namespace sgpp
