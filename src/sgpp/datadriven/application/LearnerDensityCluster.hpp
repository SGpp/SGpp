/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Florian Zipperle (florian.zipperle@tum.de)

#ifndef LEARNERDENSITYCLUSTER_HPP
#define LEARNERDENSITYCLUSTER_HPP

#include "base/datatypes/DataVector.hpp"
#include "datadriven/application/LearnerBase.hpp"
#include "datadriven/tools/TypesDatadriven.hpp"
#include "datadriven/algorithm/DensitySystemMatrix.hpp"
#include <vector>

namespace sg {

  namespace datadriven {


    class LearnerDensityCluster: public sg::datadriven::LearnerBase {
      protected:
      //DensitySystemMatrix* DMatrix_ = NULL;
        std::vector<int>* component_;
      public:

    	LearnerDensityCluster();
    	//~LearnerDensityCluster();

    	//void construct();
    	//void deserialize(const char* filename);
    	//void serialize(const char* filename);
    	//void train();
    	//void assignCluster(sg::base::DataVector point);
    	//void assignCluster(sg::base::DataMatrix points);
    	sg::datadriven::LearnerTiming train(sg::base::DataMatrix& testDataset,  sg::base::DataVector& classes,
    						  const sg::base::RegularGridConfiguration& GridConfig, const sg::solver::SLESolverConfiguration& SolverConfig,
    						  const double lamda);

    	sg::datadriven::DMSystemMatrixBase* createDMSystem(sg::base::DataMatrix& trainDataset, double lambda);

    };

  }

}

#endif /* LEARNERDENSITYCLUSTER_HPP */
