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

	double raw_data[] = {  0.07450712,0.97926035,  0.09715328,0.22845448,  0.96487699,0.96487946,  0.23688192,0.11511521,  0.92957884,0.08138401,  0.93048735,0.93014054,  0.03629434,0.71300796,  0.24126233,0.41565687,  0.34807533,0.5471371 ,  0.36379639,0.28815444,  0.71984732,0.46613355,  0.51012923,0.28628777,  0.41834259,0.51663839,  0.32735096,0.5563547 ,  0.4099042, 0.95624594,  0.40974401,0.27784173,  0.49797542,0.84134336,  0.62338174,0.81687345,  0.53132954,0.70604948,  0.30077209,0.02952919,  0.61076999,0.02570524};
	DataMatrix data(raw_data, 21,2);

	DataVector classes(2);
	classes.setAll(0);

	sg::base::RegularGridConfiguration gridConf;
	gridConf.dim_ = 2;
	gridConf.level_ = 3;
	gridConf.type_ = sg::base::Periodic;

	sg::solver::SLESolverConfiguration solvConf;
	solvConf.eps_ = 0.001;
	solvConf.maxIterations_ = 100;
	solvConf.threshold_ = -1.0;
	solvConf.type_ = sg::solver::CG;

	clust.train(data, classes, gridConf, solvConf, 0.01);
}
