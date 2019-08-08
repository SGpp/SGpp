/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Harmonica.cpp
 *
 *  Created on: Feb 2, 2018
 *      Author: Eric Koepke
 */
#include <sgpp/datadriven/datamining/modules/hpo/harmonica/Harmonica.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/solver/sle/fista/LassoFunction.hpp>
#include <sgpp/solver/sle/fista/Fista.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/harmonica/OperationMultipleEvalMatrix.hpp>

#include <vector>
#include <string>
#include <algorithm>

namespace sgpp {
namespace datadriven {

Harmonica::Harmonica(FitterFactory *fitterFactory)
    : fitterFactory(fitterFactory), configIDs(), savedScores(), configBits(), constraints() {
  fitterFactory->getConfigBits(configBits);
  freeBits = std::vector<ConfigurationBit *>(configBits);
}

void Harmonica::transformScores(const DataVector &source, DataVector &target) {
  for (size_t i = 0; i < source.size(); i++) {
    target[i] = std::pow(source[i], 0.5);
  }
}

void Harmonica::prepareConfigs(std::vector<ModelFittingBase*> &fitters,
                               int seed,
                               std::vector<std::string> &configStrings) {
  // migrate samples that fit in the new space
  size_t nOld = configIDs.size();
  // std::cout << "nOld: " << nOld << std::endl;
  size_t nAll = nOld + fitters.size();
  // std::cout << "nAll: " << nAll << std::endl;
  configIDs.resize(nAll);

  // build row of bit combinations used for the matrix and constraint reconstruction
  size_t ncols = (freeBits.size() * freeBits.size() + 5) * freeBits.size() / 6 + 1;
  // std::cout << "nBits: " << ncols << std::endl;
  // std::cout << "nConfigBits: " << configBits.size() << std::endl;
  parityrow = std::vector<std::vector<ConfigurationBit *>>(ncols);
  size_t cnt = freeBits.size();
  size_t cnt2 = (freeBits.size() + 1) * freeBits.size() / 2;
  for (size_t i = 0; i < freeBits.size(); i++) {
    parityrow[i].push_back(freeBits[i]);
    for (size_t k = i + 1; k < freeBits.size(); k++) {
      parityrow[cnt].push_back(freeBits[i]);
      parityrow[cnt].push_back(freeBits[k]);
      cnt++;
      for (size_t m = k + 1; m < freeBits.size(); m++) {
        parityrow[cnt2].push_back(freeBits[i]);
        parityrow[cnt2].push_back(freeBits[k]);
        parityrow[cnt2].push_back(freeBits[m]);
        cnt2++;
      }
    }
  }

  paritymatrix = DataMatrix(nAll, ncols);

  createRandomConfigs(freeBits.size(), configIDs, seed, nOld);

  for (size_t i = 0; i < nAll; i++) {
    setParameters(configIDs[i], i);
    if (i >= nOld) {
      fitters[i - nOld] = fitterFactory->buildFitter();
      configStrings[i - nOld] = fitterFactory->printConfig();
    }
  }
}

void Harmonica::calculateConstrainedSpace(const DataVector &transformedScores,
                                          double lambda,
                                          int shrink) {
  size_t nOld = savedScores.size();
  size_t nAll = transformedScores.size() + nOld;
  savedScores.resize(nAll);
  for (size_t i = nOld; i < nAll; i++) {
    savedScores[i] = transformedScores[i - nOld];
  }

  base::DataVector normed(savedScores);
  normed.normalize();

  // run solver
  solver::LassoFunction g{lambda};
  solver::Fista<solver::LassoFunction> fista{g};
  DataVector alpha = DataVector{paritymatrix.getNcols()};
  base::LinearGrid dummygrid(0);
  OperationMultipleEvalMatrix opMultEval{dummygrid, paritymatrix};
  fista.solve(opMultEval, alpha, normed, 100, DEFAULT_RES_THRESHOLD);

  std::vector<size_t> idx(alpha.size() - 1);
  for (size_t i = 0; i < idx.size(); i++) {
    idx[i] = i;
    // std::cout<<"Alpha: "<<i<<":"<<alpha[i]<<std::endl;   // bias term invisible
  }

  // sort indices based on comparing values in alpha
  sort(idx.begin(), idx.end(),
       [&alpha](int i1, int i2) { return fabs(alpha[i1]) > fabs(alpha[i2]); });

  size_t nBitsOld = freeBits.size();
  size_t i = 0;

  // save free Bits for moving configs to new space
  std::vector<ConfigurationBit *> freeBitsold(freeBits);

  while (freeBits.size() > nBitsOld - shrink && alpha[idx[i]] != 0) {
    int bias = -((alpha[idx[i]] > 0) - (alpha[idx[i]] < 0));
    if (addConstraint(idx[i], bias)) {
      std::cout << "Constraint added: " << parityrow[idx[i]][0]->getName();
      for (size_t k = 1; k < parityrow[idx[i]].size(); ++k) {
        std::cout << " * " << parityrow[idx[i]][k]->getName();
      }
      std::cout << " = " << bias << " (Weight: " << alpha[idx[i]] << ")" << std::endl;
    }
    i++;
  }

  std::cout << "Free bits in new space: " << freeBits.size() << std::endl;

  std::vector<int> configIDsNew{};
  DataVector newScores{};
  for (size_t k = 0; k < configIDs.size(); k++) {
    int moved = moveToNewSpace(configIDs[k], freeBitsold);
    // if in new space push back on vector
    // and save Scores
    if (moved >= 0) {
      newScores.push_back(savedScores[k]);
      configIDsNew.push_back(moved);
    }
  }
  configIDs = configIDsNew;
  savedScores = newScores;
}

void Harmonica::createRandomConfigs(size_t nBits,
                                    std::vector<int> &configIDs,
                                    int seed,
                                    size_t start) {
  std::mt19937 generator = std::mt19937(seed);
  std::uniform_int_distribution<int> distribution(0, static_cast<int>(std::pow(2, nBits) - 1));
  // std::cout << "MaxConfigs: " << std::pow(2, nBits) << std::endl;
  for (size_t i = start; i < configIDs.size(); i++) {
    configIDs[i] = distribution(generator);
    bool bUnchecked = true;
    while (bUnchecked) {
      bUnchecked = false;
      for (size_t k = 0; k < i; k++) {
        if (configIDs[i] == configIDs[k]) {
          configIDs[i] = distribution(generator);
          bUnchecked = true;
        }
      }
    }
  }
}

void Harmonica::fixConfigBits(bool resetFree) {
  if (resetFree) {
    freeBits = std::vector<ConfigurationBit *>{};
  }
  size_t nextFreeBit = 0;
  for (auto &constraint : constraints) {
    if (constraint->getOpenBits() == 2) {
      constraint->findComplex();
    }
  }
  bool resolved;
  while (nextFreeBit < configBits.size() - 1) {
    resolved = true;
    while (resolved) {
      resolved = false;
      for (auto &constraint : constraints) {
        if (constraint->getOpenBits() == 1) {
          constraint->resolve();
          resolved = true;
        }
      }
    }
    while (configBits[nextFreeBit]->getValue() != 0
        && nextFreeBit < configBits.size() - 1) {
      nextFreeBit++;
    }
    if (configBits[nextFreeBit]->getValue() == 0) {
      freeBits.push_back(configBits[nextFreeBit]);
      configBits[nextFreeBit]->setValue(1);
      if (!resetFree) {
        throw base::application_exception("Error in HPO: freeBits changed during evaluation.");
      }
    }
  }
}

void Harmonica::resetBits() {
  for (auto &bit : configBits) {
    bit->reset();
  }
  for (auto &constraint : constraints) {
    constraint->reset();
  }
}

void Harmonica::setParameters(int configID, size_t matrixrow) {
  resetBits();
  for (auto &bit : freeBits) {
    bit->setValue((configID & 1) * 2 - 1);
    configID = configID >> 1;
  }
  fixConfigBits(false);
  for (size_t i = 0; i < parityrow.size(); i++) {
    int tmp = 1;
    for (auto &bit : parityrow[i]) {
      tmp = tmp * bit->getValue();
    }
    paritymatrix.set(matrixrow, i, tmp);
  }

  fitterFactory->setHarmonica();
}

bool Harmonica::addConstraint(size_t idx, int bias) {
  constraints.push_back(std::make_unique<ConfigurationRestriction>(parityrow[idx], bias));
  for (auto &bit : parityrow[idx]) {
    bit->addConstraint(constraints.back().get());
    // std::cout << "Adding bit " << bit->getName() << " from constraint:" << idx << std::endl;
  }
  resetBits();
  fixConfigBits(true);
  if (!checkConstraints()) {
    constraints.pop_back();
    for (auto &bit : parityrow[idx]) {
      bit->removeLastConstraint();
      // std::cout << "Removing bit from constraint:" << idx << std::endl;
    }
    return false;
  }
  return true;
}

bool Harmonica::checkConstraints() {
  for (auto &constraint : constraints) {
    if (!constraint->check()) {
      return false;
    }
  }
  return true;
}

int Harmonica::moveToNewSpace(int configID, std::vector<ConfigurationBit *> oldFreeBits) {
  resetBits();
  for (auto &bit : oldFreeBits) {
    bit->setValue((configID & 1) * 2 - 1);
    configID = configID >> 1;
  }
  fixConfigBits(false);
  if (!checkConstraints()) {
    return -1;
  }
  int v = 0;
  int m = 1;
  for (auto &bit : freeBits) {
    v = v + m * ((bit->getValue() + 1) / 2);
    m = m * 2;
  }
  return v;
}
} /* namespace datadriven */
} /* namespace sgpp */
