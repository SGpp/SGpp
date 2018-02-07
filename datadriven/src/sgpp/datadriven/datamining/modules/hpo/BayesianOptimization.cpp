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



namespace sgpp {
namespace datadriven {

BayesianOptimization::BayesianOptimization(double firstvalue)
:kernelmatrix(1,1,1), kernelinv(1,1,1), transformedOutput(1), screwedvar(false), testknew(1){
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

	base::DataVector tmp(knew.size());
	//optimization::sle_solver::Eigen solver{};
    //optimization::FullSLE sle(kernelmatrix);

	//solver.solve(*sle, knew, tmp);
	optimization::sle_solver::Eigen solver{}; //100, 10e-5, tmp
	// itsolver.setTolerance(0.1);
	solver.solve(*sle, knew, tmp);
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
	return 1000000*((mean - (bestsofar-0.001))*(0.5+0.5*std::erf(z/1.41))-var*0.4*std::exp(-0.5*z*z));

}



void BayesianOptimization::updateGP(base::DataVector knew, base::DataVector y){
    optimization::sle_solver::BiCGStab solver{}; //Eigen
    base::DataVector kold(testknew);
    //kold.resize(knew.size()-1);
    base::DataVector tmp(kold.size());
    	solver.solve(*sle, kold, tmp);
    	 std::cout << kold.toString() << std::endl;
    	  base::DataVector check(kold.size());
    	  kernelmatrix.mult(tmp, check);
    	  std::cout << check.toString() << std::endl;
    	  std::cout << "Var: "<<kold.dotProduct(tmp)<<std::endl;
    	  std::cout << tmp.toString() << std::endl;

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

    bool okay = solver.solve(*sle, y, transformedOutput);
   // optimization::Printer::getInstance().disableStatusPrinting();



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

} /* namespace datadriven */
} /* namespace sgpp */
