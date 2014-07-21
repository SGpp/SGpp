/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Florian Zipperle (florian.zipperle@tum.de)

#include "base/grid/Grid.hpp"
#include "LearnerDensityCluster.hpp"
#include "solver/sle/ConjugateGradients.hpp"
#include "pde/operation/PdeOpFactory.hpp"
#include "base/operation/OperationMatrix.hpp"
#include "base/operation/OperationIdentity.hpp"
#include "base/exception/data_exception.hpp"
#include "base/datatypes/DataVector.hpp"

#include <iostream>
#include <algorithm>


#include "alglib/stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "alglib/alglibmisc.h"

using namespace alglib;

#include <boost/config.hpp>
#include <iostream>
#include <algorithm>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>

namespace sg {

  namespace datadriven {

  	LearnerDensityCluster::LearnerDensityCluster() : LearnerBase(false, false){

  	}

  	LearnerDensityCluster::~LearnerDensityCluster() {
  		delete component_;
  	}

  	const std::string currentDateTime() {
  	    time_t     now = time(0);
  	    struct tm  tstruct;
  	    char       buf[80];
  	    tstruct = *localtime(&now);
  	    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
  	    // for more information about date/time format
  	    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

  	    return buf;
  	}

    /**
     * Learning a dataset with regular sparse grids
     *
     * @param testDataset the training dataset
     * @param classes classes corresponding to the training dataset
     * @param GridConfig configuration of the regular grid
     * @param SolverConfig configuration of the SLE solver
     * @param lamda regularization parameter lambda
     */
    LearnerTiming LearnerDensityCluster::train(sg::base::DataMatrix& trainDataset,  sg::base::DataVector& classes,
					  const sg::base::RegularGridConfiguration& GridConfig, const sg::solver::SLESolverConfiguration& SolverConfig,
					  const double lambda){

		LearnerTiming result;
		result.timeComplete_ = 0.0;
		result.timeMultComplete_ = 0.0;
		result.timeMultCompute_ = 0.0;
		result.timeMultTransComplete_ = 0.0;
		result.timeMultTransCompute_ = 0.0;
		result.timeRegularization_ = 0.0;
		result.GFlop_ = 0.0;
		result.GByte_ = 0.0;

		if(GridConfig.type_ == sg::base::Periodic){
			grid_ = sg::base::Grid::createPeriodicGrid(GridConfig.dim_);
			grid_->createGridGenerator()->regular(GridConfig.level_);
		}else{
			throw base::data_exception ("Just periodic grids are supported now.");
		}


		//Equation: (A + lI)a = f
		//A = (phi_i,phi_j)_L2  - LTwoDot
		//l = lambda
		//I = identity matrix
		//a = alpha
		//f_i = 1/M * sum_j=1^M (phi_i(x_j))
		//M = number of data points (testDataset)

		std::cout << "Start " << currentDateTime() << std::endl;

    	sg::solver::SLESolver* myCG;
    	myCG = new sg::solver::ConjugateGradients(SolverConfig.maxIterations_, SolverConfig.eps_);


    	sg::base::OperationIdentity C_;
		sg::datadriven::DensitySystemMatrix DMatrix(*grid_, trainDataset, C_, lambda);

		sg::base::DataVector rhs(grid_->getStorage()->size());
		DMatrix.generateb(rhs);
		//std::cout << rhs.toString() << std::endl;

		sg::base::DataVector alpha(grid_->getStorage()->size());
		alpha.setAll(0);

		std::cout << "Matrizen erstell t" << currentDateTime() << std::endl;

		myCG->solve(DMatrix, alpha, rhs, false, true, -1.0);
		alpha_ = new sg::base::DataVector(alpha);

		std::cout << "geloest " << currentDateTime() << std::endl;
		//std::cout << alpha.toString() << std::endl;

		//alpha_ = new sg::base::DataVector(alpha);
		// calculate for each point the estimated density
		// insert the points to the kdtree from alglib
		// calculate the 5 nearest points
		// if the difference between the points is under the threshold
		// add the points to the similarity graph from boost
		// get the result

		double *points = new double[(GridConfig.dim_ + 2) * trainDataset.getNrows()];
		for(size_t i = 0; i < trainDataset.getNrows();i++){
			sg::base::DataVector point(GridConfig.dim_);
			trainDataset.getRow(i,point);
			double y = grid_->eval(alpha,point);
			points[(GridConfig.dim_+2)*(i+1)-2] = double(i);
			points[(GridConfig.dim_+2)*(i+1)-1] = y;

			std::copy(point.getPointer(), point.getPointer() + GridConfig.dim_, points+(i)*(GridConfig.dim_+2));
		}


		ae_int_t nx = GridConfig.dim_;
		ae_int_t ny = 2;
		ae_int_t normtype = 3;
		kdtree kdt;
		real_1d_array x = real_1d_array();
		real_2d_array r = "[[]]";
		real_1d_array dist = "[]";

		real_2d_array a = real_2d_array();
		a.setcontent(trainDataset.getNrows(), GridConfig.dim_ + 2, points);
		delete points;

		std::cout << "alglib array erstellt" << currentDateTime() << std::endl;

		//printf("%s\n", a.tostring(1).c_str());
		kdtreebuild(a, nx, ny, normtype, kdt);

		// query all neighbors
		typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;
		Graph G;

		double eps = 2;
		size_t k = 5;
		for(size_t i = 0; i < trainDataset.getNrows();i++){
			sg::base::DataVector point(GridConfig.dim_);
			trainDataset.getRow(i,point);

			x.setcontent(GridConfig.dim_, point.getPointer());
			kdtreequeryknn(kdt, x, k);
			kdtreequeryresultsxy(kdt, r);

			kdtreequeryresultsdistances(kdt, dist);
			// check the threshold
			//printf("----------------------------\n%d\n%d\n", int(i) , int(i));
			for(size_t j = 1; j < k; j++){
				double neighborY = r[j][GridConfig.dim_+1];
				if(neighborY >= eps || r[0][GridConfig.dim_+1] >= eps){
					//add to boost graph
					//printf("(%d) <-> (%d)\n",int(r[0][GridConfig.dim_]), int(r[j][GridConfig.dim_]));
					boost::add_edge(int(r[0][GridConfig.dim_]), int(r[j][GridConfig.dim_]), G);
				}
				//printf("%d %.5f\n",int(r[j][GridConfig.dim_]), dist[j]);
			}
			boost::add_edge(i, i, G);
			//printf("\n");
		}
		std::vector<int> component(boost::num_vertices(G));
		int num = boost::connected_components(G, &component[0]);

		std::cout << "component erstellt " << currentDateTime() << std::endl;
		/*std::vector<int>::size_type i;
		std::cout << "Total number of components: " << num << std::endl;
		for (i = 0; i != component.size(); ++i){
			std::cout << "Vertex " << i <<" is in component " << component[i] << std::endl;
		}
		std::cout << std::endl;*/

		if(component_ != NULL)
			delete component_;
		component_ = new sg::base::DataVector(component);
		//printf("done");
		delete myCG;
    	return result;
    }

    sg::datadriven::DMSystemMatrixBase* LearnerDensityCluster::createDMSystem(sg::base::DataMatrix& trainDataset, double lambda){
    	return NULL;//DMatrix_;
    }

    sg::base::DataVector LearnerDensityCluster::getComponent(){
    	return *component_;
    }
  }
}

