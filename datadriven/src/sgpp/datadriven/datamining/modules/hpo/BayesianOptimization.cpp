/*
 * BayesianOptimization.cpp
 *
 *  Created on: Feb 2, 2018
 *      Author: polarbart
 */

#include <sgpp/optimization/sle/solver/Eigen.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/BayesianOptimization.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <iostream>
#include <sgpp/optimization/sle/solver/BiCGStab.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <cmath>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/sle/solver/GaussianElimination.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/optimizer/unconstrained/MultiStart.hpp>


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
  bestsofar = rawScores.min();

  CholeskyDecomposition();
  transformedOutput = base::DataVector(rawScores);
  solveCholeskySystem(transformedOutput);
}

double BayesianOptimization::transformScore(double original) {
  return -1/(1+original);
}

double BayesianOptimization::kernel(double distance) {
  double kernelwidth = 1.5; //EDIT: metaparameter and based on dimensions?
  return exp(- distance/ 2 * kernelwidth);
}

double BayesianOptimization::acquisitionEI(base::DataVector knew, double kself, double bestsofar){
  double mean = knew.dotProduct(transformedOutput);

  base::DataVector tmp(knew); //size
  //optimization::sle_solver::Eigen solver{};
  //optimization::FullSLE sle(kernelmatrix);

  //solver.solve(*sle, knew, tmp);
  optimization::sle_solver::GaussianElimination solver{}; //100, 10e-5, tmp
  solveCholeskySystem(tmp);
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
  return ((mean - (bestsofar-0.001))*(0.5+0.5*std::erf(-z/1.41))-var*0.4*std::exp(-0.5*z*z)); //erf z

}

BOConfig* BayesianOptimization::main(BOConfig& prototype) {
  BOConfig nextconfig(prototype);
  optimization::WrapperScalarFunction wrapper(prototype.getContSize(), std::bind(acquisitionOuter, this, std::placeholders::_1));
  double min = 1.0/0; //infinity
  BOConfig bestConfig;
  //EDIT: What if no discrete exist?
  while(nextconfig.nextDisc()){
    for(auto& config: allConfigs){
      config.calcDiscDistance(nextconfig);
    }
    optimization::optimizer::MultiStart optimizer(wrapper);
    optimizer.optimize();
    if (optimizer.getOptimalValue() < min) {
      min = optimizer.getOptimalValue();
      bestConfig = BOConfig(nextconfig);
      bestConfig.setCont(optimizer.getOptimalPoint());
    }
  }
  std::cout << "Acquistion: " << min << std::endl;
  allConfigs.push_back(bestConfig);
  return &allConfigs.back();
}


double BayesianOptimization::acquisitionOuter(const base::DataVector & inp) {
  base::DataVector kernelrow(allConfigs.size());
  for (int i = 0; i < allConfigs.size(); i++) {
    kernelrow[i] = kernel(allConfigs[i].getTotalDistance(inp)); // divided by 2
    if(kernelrow[i]==1){ //EDIT: good???
      return 1.0/0;
    }
    // std::cout << "Kernel value: "<<exp((-squaresum[i] - std::pow(tmp.l2Norm(), 2)) / 2) <<std::endl;
  }
  return acquisitionEI(kernelrow, 1, bestsofar);  //EDIT: ascend or descent?, kself + noise?
}

/*
double BayesianOptimization::mean(base::DataVector knew){
	return knew.dotProduct(transformedOutput);
}

double BayesianOptimization::var(base::DataVector knew, double kself){
	base::DataVector tmp(knew.size());
	// kernelinv.mult(knew, tmp);
    optimization::FullSLE sle(kernelmatrix);
    optimization::sle_solver::Eigen solver{};
    solver.solve(sle, knew, tmp);

  if(knew.dotProduct(tmp) > 1 || knew.dotProduct(tmp)<0){
      //std::cout << "Error: wrong variance: "<<knew.dotProduct(tmp) <<", knew:"
      //         <<knew.toString()<<", temp: "<<tmp.toString()<<std::endl;
    screwedvar = true;
    return 0;
  }
	return kself - knew.dotProduct(tmp);
}

double BayesianOptimization::acquisitionPI(base::DataVector knew, double kself, double bestsofar){
    if(var(knew,kself)==0){
      return 0;
    }
    return (mean(knew) - (bestsofar-1))/var(knew, kself); //-5
}
 */




void BayesianOptimization::updateGP(){
  double noise = 1e-4; //EDIT: metaparameter
  size_t size = kernelmatrix.getNcols();
  //kernelmatrix.resize(size + 1, size + 1); //EDIT: does this work? did I have problems before? IT does not!
  kernelmatrix.appendRow();
  kernelmatrix.appendCol(base::DataVector(size+1));
  for (size_t i = 0; i < size; ++i) {
    double tmp = kernel(allConfigs[i].getDistance(allConfigs.back()));
    kernelmatrix.set(size, i, tmp);
    kernelmatrix.set(i, size, tmp);
  }
  kernelmatrix.set(size,size, 1+noise);


  rawScores.push_back(transformScore(allConfigs.back().getScore()));
  bestsofar = std::min(rawScores.back(), bestsofar);

  CholeskyDecomposition();
  transformedOutput = base::DataVector(rawScores);
  solveCholeskySystem(transformedOutput);


  // optimization::sle_solver::GaussianElimination solver{}; //Eigen
  // sle = std::make_unique<optimization::FullSLE>(kernelmatrix);
  // bool okay = solver.solve(*sle, y, transformedOutput);
  // std::cout << "Solver okay: " << okay << std::endl;


  std::cout << "Maxofmax: "<<maxofmax<<std::endl;
  maxofmax = 0;
  base::DataVector check2(transformedOutput.size());
  kernelmatrix.mult(transformedOutput, check2);
  check2.sub(rawScores);
  double max = check2.maxNorm();
  std::cout << "Max: "<<max<<std::endl;
   //std::cout<<transformedOutput.toString()<<std::endl;

  std::cout << "Var Screwed: " << screwedvar << std::endl;
  screwedvar = false;

}

void BayesianOptimization::CholeskyDecomposition(){
  size_t n = kernelmatrix.getNrows();
  gleft = base::DataMatrix(n, n, 0);

  for(size_t i=0; i < n; i++){
  for(size_t j=0; j <= i; j++){
  double sum = kernelmatrix.get(i,j);
  for(size_t k=0; k < j; k++){
    sum = sum - gleft.get(i, k) * gleft.get(j, k);
  }
  if(i > j){
    double value = sum/gleft.get(j,j);
    gleft.set(i, j, value);
    //gright.set(j, i, value);
  }else if(sum > 0){
    double value = std::sqrt(sum);
    gleft.set(i, i, value);
    //gright.set(i, i, value);
  }else {
    std::cout << "Error: CholDecom failed." << std::endl;
    gleft.set(i,i,10e-8); //EDIT: experimental
  }
  }
  }

}

void BayesianOptimization::solveCholeskySystem(base::DataVector& x){
  // base::DataVector x(b);
  // std::cout<<"Test Point 1"<<std::endl;
  for(size_t i = 0; i < x.size(); i++){
    x[i] = x[i] / gleft.get(i,i);
    for(size_t k = i+1; k < x.size();k++){
      x[k] = x[k] - gleft.get(k, i)*x[i];
    }
  }
  //std::cout<<"Test Point 2"<<std::endl;

  for(int i = x.size()-1; i >= 0; i--){
    x[i] = x[i] / gleft.get(i,i);
    for(int k = i-1; k >= 0;k--){
      x[k] = x[k] - gleft.get(i, k)*x[i];
    }
  }
  // std::cout<<"Test Point 3"<<std::endl;
}




} /* namespace datadriven */
} /* namespace sgpp */
