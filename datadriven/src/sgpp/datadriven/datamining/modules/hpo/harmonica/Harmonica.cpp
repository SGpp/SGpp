/*
 * Harmonica.cpp
 *
 *  Created on: Feb 2, 2018
 *      Author: Eric Koepke
 */

#include <sgpp/optimization/sle/solver/Eigen.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/harmonica/Harmonica.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <iostream>
#include <sgpp/optimization/sle/solver/BiCGStab.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <cmath>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/sle/solver/GaussianElimination.hpp>
#include <random>
#include <sgpp/solver/sle/fista/LassoFunction.hpp>
#include <sgpp/solver/sle/fista/Fista.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include "OperationMultipleEvalMatrix.hpp"


namespace sgpp {
namespace datadriven {

Harmonica::Harmonica(FitterFactory* fitterFactory)
:fitterFactory(fitterFactory), configIDs(), savedScores(), configBits(), constraints() {
  fitterFactory->getConfigBits(configBits);
  freeBits = std::vector<ConfigurationBit *>(configBits);
  for (auto &bit : freeBits) {
    std::cout << "freeBits: " << bit->getName() << std::endl;
  }
  for (auto &bit : configBits) {
    std::cout << "configBits: " << bit->getName() << std::endl;
  }
}

/*
void Harmonica::transformScores(const DataVector& source, DataVector& target){
  std::vector<int> idx(source.size());
  for(int i=0; i<idx.size(); i++){
    idx[i] = i;
  }
  // sort indices based on comparing values in source
  sort(idx.begin(), idx.end(),
       [&source](int i1, int i2) {return source[i1] < source[i2];});
  for(int i=0; i<idx.size(); i++){
    target[idx[i]] = i;
    if(i>50){
      target[idx[i]] = 50;
    }
  }
}
*/

void Harmonica::transformScores(const DataVector& source, DataVector& target){
  for(int i=0;i<source.size();i++){
    target[i] = std::pow(source[i],0.5);
  }
}

std::vector<int> * Harmonica::prepareConfigs(std::vector<std::unique_ptr<ModelFittingBase>> &fitters, int seed,
                                             std::vector<std::string> &configStrings) {

  //migrate samples that fit in the new space
  size_t nOld = configIDs.size();
  std::cout<<"nOld: "<<nOld << std::endl;
  size_t nAll = nOld + fitters.size();
  std::cout<<"nAll: "<<nAll << std::endl;
  configIDs.resize(nAll);

  // get configured models for n samples (call fitterfactory)
  // EDIT: deleted freeBit assignment here, add in constructor
  size_t ncols = (freeBits.size()*freeBits.size()+5)*freeBits.size()/6 +1;
  std::cout<<"nBits: "<<ncols<<std::endl;
  std::cout<<"nConfigBits: "<<configBits.size()<<std::endl;
  parityrow = std::vector<std::vector<ConfigurationBit*>>(ncols);
  /*for(int i = 0; i<ncols;i++){
    parityrow.emplace_back();
  }*/
  size_t cnt = freeBits.size();
  size_t cnt2 =(freeBits.size() +1)*freeBits.size()/2;
  for(int i = 0; i < freeBits.size(); i++){
    parityrow[i].push_back(freeBits[i]);
    for(int k = i+1; k < freeBits.size(); k++){
      parityrow[cnt].push_back(freeBits[i]);
      parityrow[cnt].push_back(freeBits[k]);
      cnt++;
      for(int m = k+1; m < freeBits.size(); m++){
        parityrow[cnt2].push_back(freeBits[i]);
        parityrow[cnt2].push_back(freeBits[k]);
        parityrow[cnt2].push_back(freeBits[m]);
        cnt2++;
      }
    }
  }

  //nAll = 1000;
  //build matrix
  // EDIT: add bias term? already in there!!!
  paritymatrix = DataMatrix(nAll, ncols);

  createRandomConfigs(freeBits.size(), configIDs, 42, nOld); //EDIT: seed
  /*std::vector<int> shuffleConfigs(configIDs);
  std::mt19937 sgenerator(seed);
  std::shuffle(shuffleConfigs.begin(), shuffleConfigs.end(), sgenerator);
  configIDs.resize(nAll);
  for (int j = 0; j < nAll; ++j) {
    configIDs[j] = shuffleConfigs[j];
  }*/

  std::cout << "ConfigIDs: " << configIDs[0] << ", " << configIDs[1] << std::endl;

  for (int i = 0; i < nAll; i++) {
    setParameters(configIDs[i], i);
    if(i>=nOld) {
      fitters[i-nOld].reset(fitterFactory->buildFitter());
      configStrings[i-nOld] = fitterFactory->printConfig();
    }
  }
  return &configIDs;
}


//EDIT: handle zero bias constraints

void Harmonica::calculateConstrainedSpace(const DataVector& transformedScores, double lambda, int shrink){
  size_t nOld = savedScores.size();
  size_t nAll = transformedScores.size()+nOld;
  savedScores.resize(nAll); //EDIT: is this working?
  for(size_t i = nOld; i<nAll;i++){
    savedScores[i] = transformedScores[i-nOld];
  }

  // run solver
  solver::LassoFunction g{lambda};
  solver::Fista<solver::LassoFunction> fista{g};
  DataVector alpha = DataVector{paritymatrix.getNcols()};
  base::LinearGrid dummygrid(0);
  OperationMultipleEvalMatrix opMultEval{dummygrid, paritymatrix}; //EDIT: This okay?
  fista.solve(opMultEval, alpha, savedScores, 100, DEFAULT_RES_THRESHOLD);

  std::vector<int> idx(alpha.size()-1);
  for(int i=0; i<idx.size(); i++){
    idx[i] = i;
    // std::cout<<"Alpha: "<<i<<":"<<alpha[i]<<std::endl; //bias term invisible
  }

  // sort indices based on comparing values in alpha
  sort(idx.begin(), idx.end(),
       [&alpha](int i1, int i2) {return fabs(alpha[i1]) > fabs(alpha[i2]);});

  size_t nBitsOld = freeBits.size();
  int i = 0;

  //save free Bits for moving configs to new space
  std::vector<ConfigurationBit*> freeBitsold(freeBits);

  while(freeBits.size() > nBitsOld-shrink && alpha[idx[i]] != 0){
    addConstraint(idx[i], -((alpha[idx[i]]>0)-(alpha[idx[i]]<0))); // -lround(alpha[idx[i]]/fabs(alpha[idx[i]]))
    std::cout<<"Constraint number:"<<idx[i]<<","<<-((alpha[idx[i]]>0)-(alpha[idx[i]]<0))<<","<<alpha[idx[i]]<<", freebits: "<<freeBits.size()<<std::endl;
    i++;
  }

  std::vector<int> configIDsNew{};
  DataVector newScores{};
  for(int i=0; i< configIDs.size();i++){
    int moved = moveToNewSpace(configIDs[i], freeBitsold);
    // if in new space push back on vector
    // and save Scores
    if(moved >= 0){
      newScores.push_back(savedScores[i]);
      configIDsNew.push_back(moved);
    }
  }
  configIDs = configIDsNew; //EDIT: is this working?
  savedScores = newScores;
}



void Harmonica::createRandomConfigs(size_t nBits, std::vector<int>& configIDs, int seed, size_t start) {
  std::mt19937 generator = std::mt19937(seed);
  std::uniform_int_distribution<int> distribution(0, static_cast<int>(std::pow(2, nBits) - 1));
  std::cout << "MaxConfigs: " << std::pow(2, nBits) << std::endl;
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

bool Harmonica::fixConfigBits() { //EDIT: einfacherer Approach, nur über constraints gehen
  //EDIT: reset bits and constraints
  //std::vector<ConfigurationBit*> freeBitsn{};
  int nextFreeBit = 0;
  bool changedFreeBits = false;
  bool resolved = true;
  while(nextFreeBit < configBits.size()-1) {
    resolved = true;
    while (resolved) {
      resolved = false;
      for (auto &constraint: constraints) {
        if (constraint->getOpenBits() == 1) {
          constraint->resolve();
          resolved = true;
        }
      }
    }
    while (configBits[nextFreeBit]->getValue() != 0 && nextFreeBit < configBits.size() - 1) { //EDIT: endpoint
      nextFreeBit++;
    }
    if (configBits[nextFreeBit]->getValue() == 0) {
      freeBits.push_back(configBits[nextFreeBit]);
      configBits[nextFreeBit]->setValue(1);
      changedFreeBits = true;
    }
  }
  return changedFreeBits;
}

void Harmonica::resetBits() {
  for (auto &bit : configBits) {
    bit->reset();
  }
  for (auto &constraint: constraints) {
    constraint->reset();
  }
}

void Harmonica::setParameters(int configID, int matrixrow) {
  resetBits();
  for(auto& bit : freeBits){
    bit->setValue((configID&1)*2-1);
    configID = configID >> 1;
  }
  bool changed = fixConfigBits();
  if(changed){
    std::cout<<"Error: freeBits changed in setParameters."<<std::endl; //EDIT: throw exception
  }
  for(int i=0;i<parityrow.size();i++){
    int tmp = 1;
    for(auto& bit : parityrow[i]){
      tmp = tmp * bit->getValue();
    }
    paritymatrix.set(matrixrow, i, tmp);
  }

  fitterFactory->setHarmonica();

  //EDIT: set Parameters in fitter and parameter classes
}

void Harmonica::addConstraint(int idx, int bias){
  constraints.push_back(std::make_unique<ConfigurationRestriction>(parityrow[idx], bias));
  for(auto &bit : parityrow[idx]){
    bit->addConstraint(constraints.back().get());
    std::cout<<"Adding bit "<<bit->getName()<<" from constraint:"<<idx<<std::endl;
  }
  freeBits = std::vector<ConfigurationBit*>{};
  resetBits();
  fixConfigBits();
  if(!checkConstraints()){
     //EDIT: revert constraint
     constraints.pop_back();
     for(auto &bit : parityrow[idx]){
       bit->removeLastConstraint();
       std::cout<<"Removing bit from constraint:"<<idx<<std::endl;
     }
  }
}

bool Harmonica::checkConstraints() {
  for (auto &constraint: constraints) {
    if(!constraint->check()){
      return false;
    }
  }
  return true;
}

int Harmonica::moveToNewSpace(int configID, std::vector<ConfigurationBit*> oldFreeBits) {
  resetBits();
  for(auto& bit : oldFreeBits){
    bit->setValue((configID&1)*2-1);
    configID = configID >> 1;
  }
  bool changed = fixConfigBits();
  if(changed){
    std::cout<<"Error: freeBits changed in moveToNewSpace."<<std::endl; //EDIT: throw exception
  }
  if(!checkConstraints()){
    return -1;
  }
  int v = 0;
  int m = 1;
  for(auto& bit : freeBits){
    v = v+m*((bit->getValue()+1)/2);
    m = m*2;
  }
  return v;
}

} /* namespace datadriven */
} /* namespace sgpp */
