/*
 * Harmonica.cpp
 *
 *  Created on: Feb 2, 2018
 *      Author: polarbart
 */

#include <sgpp/optimization/sle/solver/Eigen.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/Harmonica.hpp>
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
:fitterFactory(fitterFactory), nBits(0), paritymatrix(0,0), configIDs(), savedScores(){}

void Harmonica::transformScores(const DataVector& source, DataVector& target){
  for(int i=0;i<source.size();i++){
    target[i] = std::pow(source[i],0.5);
  }
}

void Harmonica::prepareConfigs(std::vector<ModelFittingBase*>& fitters) {

  //migrate samples that fit in the new space
  std::vector<int> configIDsNew{};
  DataVector newScores{};
  for(int i=0; i< configIDs.size();i++){
    int moved = fitterFactory->moveToNewSpace(configIDs[i]);
    // if in new space push back on vector
    // and save Scores
    if(moved >= 0){
      newScores.push_back(savedScores[i]);
      configIDsNew.push_back(configIDs[i]);
    }
  }
  size_t nOld = configIDsNew.size();
  size_t nAll = nOld + fitters.size();
  configIDsNew.resize(nAll);
  configIDs = configIDsNew; //EDIT: is this working?

  // get configured models for n samples (call fitterfactory)
  nBits = fitterFactory->buildParity();

  //build matrix
  // EDIT: add bias term? already in there!!!
  paritymatrix = DataMatrix(nAll, (nBits * nBits + 5) * nBits / 6 + 1);

  createRandomConfigs(nBits, configIDs, 42, nOld); //EDIT: seed

  for (int i = 0; i < nAll; i++) {
    fitterFactory->setHarmonica(configIDs[i], i, paritymatrix);
    if(i>=nOld) {
      fitters[i-nOld] = fitterFactory->buildFitter();
    }
  }

}

//EDIT: handle zero bias constraints

void Harmonica::calculateConstrainedSpace(const DataVector& transformedScores, int lambda=2, int shrink=4){
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
  OperationMultipleEvalMatrix opMultEval{*new base::LinearGrid(0), paritymatrix}; //EDIT: grid l�schen
  fista.solve(opMultEval, alpha, savedScores, 100, DEFAULT_RES_THRESHOLD);

  std::vector<int> idx(alpha.size()-1);
  for(int i=0; i<idx.size(); i++){
    idx[i] = i;
    std::cout<<"Alpha: "<<i<<":"<<alpha[i]<<std::endl; //bias term invisible
  }

  // sort indices based on comparing values in alpha
  sort(idx.begin(), idx.end(),
       [&alpha](int i1, int i2) {return fabs(alpha[i1]) > fabs(alpha[i2]);});

  int freebits = nBits;
  int i = 0;
  //for(int i = 0; i < shrink; i++){
  while(freebits > nBits-shrink){
    freebits = fitterFactory->addConstraint(idx[i], -lround(alpha[idx[i]]/fabs(alpha[idx[i]]))); //(alpha[idx[i]]>0)-(alpha[idx[i]]<0)
    std::cout<<"Constraint number:"<<idx[i]<<","<<-lround(alpha[idx[i]]/fabs(alpha[idx[i]]))<<","<<alpha[idx[i]]<<std::endl;
    i++;
  }


}

void Harmonica::createRandomConfigs(int nBits, std::vector<int>& configIDs, int seed, int start) {
  std::mt19937 generator = std::mt19937(seed);
  std::uniform_int_distribution<int> distribution(0, std::pow(2, nBits) - 1);
  std::cout << "MaxConfigs: " << std::pow(2, nBits) << std::endl;
  for (int i = start; i < configIDs.size(); i++) {
    configIDs[i] = distribution(generator);
    bool bUnchecked = true;
    while (bUnchecked) {
      bUnchecked = false;
      for (int k = 0; k < i; k++) {
        if (configIDs[i] == configIDs[k]) {
          configIDs[i] = distribution(generator);
          bUnchecked = true;
        }
      }
    }
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
