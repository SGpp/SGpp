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


namespace sgpp {
namespace datadriven {

BayesianOptimization::BayesianOptimization(double firstvalue)
:kernelmatrix(1,1,1), kernelinv(1,1,1), transformedOutput(1), screwedvar(false),
 testknew(1), maxofmax(0), gleft(1,1,1){
	transformedOutput[0]=firstvalue;
	testknew[0] = 1;
    sle = std::make_unique<optimization::FullSLE>(kernelmatrix);
}

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
		if(!screwedvar){
			testknew = base::DataVector(knew);
		}
	    screwedvar = true;
	    return 0;
	}

	double z = (mean - (bestsofar-0.001))/var;
	return ((mean - (bestsofar-0.001))*(0.5+0.5*std::erf(-z/1.41))-var*0.4*std::exp(-0.5*z*z)); //erf z

}



void BayesianOptimization::updateGP(base::DataVector knew, base::DataVector y){
    optimization::sle_solver::GaussianElimination solver{}; //Eigen
   /* base::DataVector kold(testknew);
    //kold.resize(knew.size()-1);
    base::DataVector tmp(kold.size());
  base::DataVector tmp2(kold);
  solveCholeskySystem(kold, tmp2);

  solver.solve(*sle, kold, tmp);
    	 std::cout << "kold: "<<kold.toString() << std::endl;
    	  base::DataVector check(kold.size());
    	  kernelmatrix.mult(tmp2, check);
    	  std::cout <<"Reverse mult: "<< check.toString() << std::endl;
    	  std::cout << "Var: "<<kold.dotProduct(tmp2)<<std::endl;
  std::cout <<"Cholesky reverse: "<< tmp2.toString() << std::endl;

  std::cout <<"Correct reverse: "<< tmp.toString() << std::endl;
    std::cout <<"Gram Matrix: "<< kernelmatrix.toString() << std::endl;
  base::DataMatrix gprod(kernelmatrix.getNrows(),kernelmatrix.getNrows(),0);
  for(int i=0;i<gprod.getNrows();i++){
    base::DataVector tmp(gprod.getNrows());
    gleft.getRow(i, tmp);
    base::DataVector res(gprod.getNrows());
    gleft.mult(tmp, res);
    gprod.setColumn(i, res);
  }
  std::cout <<"G*G^T: "<< gprod.toString() << std::endl;

*/

  std::cout << "Maxofmax: "<<maxofmax<<std::endl;
  maxofmax = 0;

  testknew = base::DataVector(knew);
	size_t idx = kernelmatrix.appendRow();
	kernelmatrix.appendCol(knew);
	kernelmatrix.setRow(idx, knew);
    // kernelinv.resize(idx+1,idx+1);
  transformedOutput.resize(y.size());
    sle = std::make_unique<optimization::FullSLE>(kernelmatrix);
    // base::DataMatrix identity(idx+1, idx+1, 0);
    // for(size_t i=0; i<idx+1;i++){
    //    identity.set(i,i,1);
    // }
  // std::cout << "vor inv solve" << std::endl;
   // bool okay = solver.solve(sle, identity, kernelinv);
   // optimization::Printer::getInstance().enableStatusPrinting();

  CholeskyDecomposition();
  transformedOutput = base::DataVector(y);
  solveCholeskySystem(transformedOutput);

    // bool okay = solver.solve(*sle, y, transformedOutput);


   // optimization::Printer::getInstance().disableStatusPrinting();

  base::DataVector check2(transformedOutput.size());
  kernelmatrix.mult(transformedOutput, check2);
  check2.sub(y);
  double max = check2.maxNorm();
  std::cout << "Max: "<<max<<std::endl;
   //std::cout<<transformedOutput.toString()<<std::endl;


  // std::cout << "Solver okay: " << okay << std::endl;
  std::cout << "Var Screwed: " << screwedvar << std::endl;
  screwedvar = false;


   /* if(!okay){
    std::cout << knew.toString() << std::endl;

    kernelmatrix.toFile("kernelmatrix.txt");

   std::cout << kernelmatrix.toString() << std::endl;
   throw base::data_exception{"Kernelmatrix invalid."};
  } */
  // std::cout  << std::endl;

  // std::cout << kernelinv.toString() << std::endl;

	// kernelinv.mult(y, transformedOutput);
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
