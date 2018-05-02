/*
 * BayesianOptimization.cpp
 *
 *  Created on: Feb 2, 2018
 *      Author: Eric Koepke
 */

#include <sgpp/optimization/sle/solver/Eigen.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/bo/BayesianOptimization.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <iostream>
#include <sgpp/optimization/sle/solver/BiCGStab.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <cmath>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/sle/solver/GaussianElimination.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/optimizer/unconstrained/MultiStart.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>


namespace sgpp {
namespace datadriven {

BayesianOptimization::BayesianOptimization(const std::vector<BOConfig>& initialConfigs)
:kernelmatrix(initialConfigs.size(),initialConfigs.size()), gleft(), transformedOutput(), rawScores(initialConfigs.size()),
 screwedvar(false), maxofmax(0), allConfigs(initialConfigs) {
  double noise = 1e-4; //: metaparameter
  for (size_t i = 0; i < allConfigs.size(); ++i) {
    for (size_t k = 0; k < i; ++k) {
      double tmp = kernel(allConfigs[i].getDistance(allConfigs[k]));
      kernelmatrix.set(k, i, tmp);
      kernelmatrix.set(i, k, tmp);
    }
    kernelmatrix.set(i, i, 1+noise);
    rawScores[i] = transformScore(allConfigs[i].getScore());
  }
 // bestsofar = rawScores.min();
  base::DataVector normed(rawScores);
  normed.normalize();
  normed.sub(base::DataVector(normed.size(), normed.sum()/normed.size()));
  bestsofar = normed.min();
  //transformedOutput = base::DataVector(meaned);
  transformedOutput = base::DataVector(normed);

  decomposeCholesky(kernelmatrix, gleft);
  //transformedOutput = base::DataVector(rawScores);
  solveCholeskySystem(gleft, transformedOutput);
  scales = base::DataVector(3,1); //EDIT: dimensions
}

double BayesianOptimization::transformScore(double original) {
  return -1/(1+original);
}

double BayesianOptimization::kernel(double distance) {
  double kernelwidth = 8.0; //EDIT: metaparameter and based on dimensions? 1.5
  return exp(- distance * kernelwidth / 2 );
}

double BayesianOptimization::acquisitionEI(base::DataVector knew, double kself, double bestsofar){
  double mean = knew.dotProduct(transformedOutput);

  base::DataVector tmp(knew); //size
  //optimization::sle_solver::Eigen solver{};
  //optimization::FullSLE sle(kernelmatrix);

  //solver.solve(*sle, knew, tmp);
  optimization::sle_solver::GaussianElimination solver{}; //100, 10e-5, tmp
  solveCholeskySystem(gleft, tmp);
  //solver.solve(*sle, knew, tmp);
  base::DataVector check(tmp.size());
  kernelmatrix.mult(tmp, check);
  check.sub(knew);
  double max = check.maxNorm();
  maxofmax = std::fmax(max, maxofmax);
  double var = kself - knew.dotProduct(tmp);
  if(var > 1 || var<0){
    //std::cout << "Error: wrong variance: "<<var;
    // <<", knew:"<<knew.toString()<<", temp: "<<tmp.toString()<<std::endl;
    screwedvar = true;
    return 0;
  }

  double z = (mean - (bestsofar-0.001))/var;
  return ((mean - (bestsofar-0.001))*(0.5+0.5*std::erf(-z/1.41))-var*0.4*std::exp(-0.5*z*z)); //EDIT: is this calculated properly?
  //return mean; //EDIT
}

BOConfig* BayesianOptimization::main(BOConfig& prototype) {
  BOConfig nextconfig(prototype);
  optimization::WrapperScalarFunction wrapper(prototype.getContSize(), std::bind(&BayesianOptimization::acquisitionOuter, this, std::placeholders::_1));
  double min = 1.0/0; //infinity
  BOConfig bestConfig;
  //EDIT: What if no discrete exist?
  do{
    for(auto& config: allConfigs){
      config.calcDiscDistance(nextconfig, scales);
    }
    optimization::optimizer::MultiStart optimizer(wrapper);
    optimizer.optimize();
    if (optimizer.getOptimalValue() < min) {
      min = optimizer.getOptimalValue();
      bestConfig = BOConfig(nextconfig);
      bestConfig.setCont(optimizer.getOptimalPoint());
    }
    /*for (int i = 0; i < 101; ++i) {
      std::ofstream myfile("C:/Users/Eric/Documents/acquisition.txt", std::ios_base::app);
      if (myfile.is_open()) {

        //myfile << "threshold,lambda,nopoints,level,basis" << std::endl;
        myfile << nextconfig.getDisc(0) << ", " << i/100.0 << ", " << acquisitionOuter(base::DataVector(1, i/100.0)) << std::endl;
      }
      myfile.close();
    }*/
  }while(nextconfig.nextDisc());
  std::cout << "Acquistion: " << min << std::endl;
  allConfigs.push_back(bestConfig);
  return &allConfigs.back();
}


double BayesianOptimization::acquisitionOuter(const base::DataVector & inp) {
  base::DataVector kernelrow(allConfigs.size());
  for (int i = 0; i < allConfigs.size(); i++) {
    kernelrow[i] = kernel(allConfigs[i].getTotalDistance(inp, scales)); // divided by 2
    if(kernelrow[i]==1){ //EDIT: good???
      return 1.0/0;
    }
    // std::cout << "Kernel value: "<<exp((-squaresum[i] - std::pow(tmp.l2Norm(), 2)) / 2) <<std::endl;
  }
  return acquisitionEI(kernelrow, 1, bestsofar);  //EDIT: ascend or descent?, kself + noise?
}

void BayesianOptimization::fitScales() {
  optimization::WrapperScalarFunction wrapper(3, std::bind(&BayesianOptimization::likelihood, this, std::placeholders::_1)); //EDIT: number of hyperparameters
  optimization::optimizer::MultiStart optimizer(wrapper, 2000, 200); //10000
  optimizer.optimize();
  std::cout<<optimizer.getOptimalPoint().toString()<<std::endl;
  std::cout<<optimizer.getOptimalValue()<<std::endl;
  base::DataVector tmp(optimizer.getOptimalPoint());
  tmp.mult(0.1); //0.3
  scales.mult(0.9); //0.7
  scales.add(tmp);
  std::cout<<scales.toString()<<std::endl;
  double noise = pow(10, -scales.back()*10); //EDIT: metaparameter 1e-4
  for (size_t i = 0; i < allConfigs.size(); ++i) {
    for (size_t k = 0; k < i; ++k) {
      double tmp = kernel(allConfigs[i].getScaledDistance(allConfigs[k], scales));
      kernelmatrix.set(k, i, tmp);
      kernelmatrix.set(i, k, tmp);
    }
    kernelmatrix.set(i, i, 1+noise);
  }
  decomposeCholesky(kernelmatrix, gleft);
}

double BayesianOptimization::likelihood(const base::DataVector &inp) {
  double noise = pow(10, -inp.back()*10);; //: metaparameter
  base::DataMatrix km(allConfigs.size(),allConfigs.size());
  for (size_t i = 0; i < allConfigs.size(); ++i) {
    for (size_t k = 0; k < i; ++k) {
      double tmp = kernel(allConfigs[i].getScaledDistance(allConfigs[k], inp));
      km.set(k, i, tmp);
      km.set(i, k, tmp);
    }
    km.set(i, i, 1+noise);
  }
  base::DataMatrix gnew;
  decomposeCholesky(km, gnew);

  base::DataVector normed(rawScores);
  normed.normalize();
  normed.sub(base::DataVector(normed.size(), normed.sum()/normed.size()));
  //transformedOutput = base::DataVector(meaned);
  base::DataVector transformed(normed);
  solveCholeskySystem(gnew, transformed);
  double tmp = 0;
  for (size_t i = 0; i < allConfigs.size(); ++i) {
    tmp += std::log(gnew.get(i,i));
  }
  return 2*tmp+normed.dotProduct(transformed);
}

void BayesianOptimization::updateGP(){
  double noise = pow(10, -scales.back()*10); //EDIT: metaparameter
  size_t size = kernelmatrix.getNcols();
  //kernelmatrix.resize(size + 1, size + 1); //EDIT: does this work? did I have problems before? IT does not!
  kernelmatrix.appendRow();
  kernelmatrix.appendCol(base::DataVector(size+1));
  for (size_t i = 0; i < size; ++i) {
    double tmp = kernel(allConfigs[i].getScaledDistance(allConfigs.back(), scales));
    kernelmatrix.set(size, i, tmp);
    kernelmatrix.set(i, size, tmp);
  }
  kernelmatrix.set(size,size, 1+noise);


  rawScores.push_back(transformScore(allConfigs.back().getScore()));
 // bestsofar = std::min(rawScores.back(), bestsofar);

  decomposeCholesky(kernelmatrix, gleft);
  base::DataVector normed(rawScores);
  normed.normalize();
  normed.sub(base::DataVector(normed.size(), normed.sum()/normed.size()));
  bestsofar = normed.min();
  //transformedOutput = base::DataVector(meaned);
  transformedOutput = base::DataVector(normed);
//  transformedOutput = base::DataVector(rawScores);
  solveCholeskySystem(gleft, transformedOutput);


  // optimization::sle_solver::GaussianElimination solver{}; //Eigen
  // sle = std::make_unique<optimization::FullSLE>(kernelmatrix);
  // bool okay = solver.solve(*sle, y, transformedOutput);
  // std::cout << "Solver okay: " << okay << std::endl;


  std::cout << "Maxofmax: "<<maxofmax<<std::endl;
  maxofmax = 0;
  base::DataVector check2(transformedOutput.size());
  kernelmatrix.mult(transformedOutput, check2);
  check2.sub(normed);
  double max = check2.maxNorm();
  std::cout << "Max: "<<max<<std::endl;
   //std::cout<<transformedOutput.toString()<<std::endl;

  std::cout << "Var Screwed: " << screwedvar << std::endl;
  screwedvar = false;

}

void BayesianOptimization::decomposeCholesky(base::DataMatrix &km, base::DataMatrix &gnew) {
  size_t n = km.getNrows();
  gnew = base::DataMatrix(n, n, 0); //EDIT: working?

  for(size_t i=0; i < n; i++){
  for(size_t j=0; j <= i; j++){
  double sum = km.get(i,j);
  for(size_t k=0; k < j; k++){
    sum = sum - gnew.get(i, k) * gnew.get(j, k);
  }
  if(i > j){
    double value = sum/gnew.get(j,j);
    gnew.set(i, j, value);
    //gright.set(j, i, value);
  }else if(sum > 0){
    double value = std::sqrt(sum);
    gnew.set(i, i, value);
    //gright.set(i, i, value);
  }else {
    std::cout << "Error: CholDecom failed." << std::endl;
    // std::cout << kernelmatrix.toString() << std::endl;
    gnew.set(i,i,10e-8); //EDIT: experimental
  }
  }
  }

}

void BayesianOptimization::solveCholeskySystem(base::DataMatrix &gmatrix, base::DataVector &x) {
  // base::DataVector x(b);
  // std::cout<<"Test Point 1"<<std::endl;
  for(size_t i = 0; i < x.size(); i++){
    x[i] = x[i] / gmatrix.get(i,i);
    for(size_t k = i+1; k < x.size();k++){
      x[k] = x[k] - gmatrix.get(k, i)*x[i];
    }
  }
  //std::cout<<"Test Point 2"<<std::endl;

  for(int i = x.size()-1; i >= 0; i--){
    x[i] = x[i] / gmatrix.get(i,i);
    for(int k = i-1; k >= 0;k--){
      x[k] = x[k] - gmatrix.get(i, k)*x[i];
    }
  }
  // std::cout<<"Test Point 3"<<std::endl;
}




} /* namespace datadriven */
} /* namespace sgpp */
