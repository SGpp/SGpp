#include <iostream>
// All SG++ headers
//#include "sgpp_base.hpp"

// Or, better!, include only those that are required
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/Grid.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/GridGenerator.hpp"
#include "base/operation/OperationEval.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "datadriven/application/LearnerDensityCluster.hpp"

using namespace std;
using namespace sg;
using namespace sg::base;



int main() {

	sg::datadriven::LearnerDensityCluster clust;

	double raw_data[] = {0.1 , 0.1,  0.2 , 0.1,  0.3 , 0.2,  0.2 , 0.3,  0.3 , 0.3,  0.7 , 0.6,  0.8 , 0.6,  0.7 , 0.7,  0.8 , 0.7 ,  0.7 , 0.8};
	DataMatrix data(raw_data, 10,2);

	DataVector classes(2);
	classes.setAll(0);

	sg::base::RegularGridConfiguration gridConf;
	gridConf.dim_ = 2;
	gridConf.level_ = 2;
	gridConf.type_ = sg::base::Periodic;

	sg::solver::SLESolverConfiguration solvConf;
	solvConf.eps_ = 0.001;
	solvConf.maxIterations_ = 1000;
	solvConf.threshold_ = -1.0;
	solvConf.type_ = sg::solver::CG;

	clust.train(data, classes, gridConf, solvConf, 0.01);
}
