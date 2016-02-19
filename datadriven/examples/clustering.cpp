#include <chrono>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/Configuration.hpp>
#include <iostream>
#include <string>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OperationDensityOCLMultiPlatform.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/KernelMult.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/OpFactory.hpp>

#include "sgpp/datadriven/tools/ARFFTools.hpp"
using namespace SGPP::base;
int main()
{

  std::string filename = "simple_test.arff";

  std::cout << "Loading file: " << filename << std::endl;
  SGPP::datadriven::Dataset data =
      SGPP::datadriven::ARFFTools::readARFF(filename);
  SGPP::base::DataMatrix& dataset = data.getData();
  std::cout<<"Loaded "<<dataset.getNcols()<<" dimensional dataset with "<<dataset.getNrows()<<" datapoints."<<std::endl;

  unsigned int dimension=dataset.getNcols(),tiefe,k;
	double lambda,treshold;
	std::cout<<"Size of Grid (3-18): ";
	std::cin>>tiefe;
	std::cout<<"Lambda (controlls smoothness of the density function. 0.01 - 0.0001 recommended.): ";
	std::cin>>lambda;
	std::cout<<"k (Number of neighbors for each datapoint. 4 - 12 recommended.): ";
	std::cin>>k;
	std::cout<<"Treshold (nodes and edges will be removed if their density is below this treshold. 0 - 0.2 recommended.): ";
	std::cin>>treshold;
	//Create Grid
	Grid* grid = Grid::createLinearGrid(dimension);
	GridGenerator* gridGen = grid->createGridGenerator();
	gridGen->regular(tiefe);
	size_t gridsize = grid->getStorage()->size();
	std::cerr<<"Grid created! Number of grid points:	 "<<gridsize<<std::endl;

	SGPP::base::DataVector alpha(gridsize);
	SGPP::base::DataVector result(gridsize);
	for(size_t i = 0; i < gridsize; i++)
	{
		alpha[i] = 1.0;
		result[i] = 0.0;
	}

	sg::solver::ConjugateGradients *solver=new sg::solver::ConjugateGradients(17,0.0001);
	SGPP::datadriven::StreamingOCLMultiPlatform::OperationDensityOCL* operation_mult=
		SGPP::datadriven::createDensityOCLMultiPlatformConfigured(*grid, dimension, lambda, "MyOCLConf.cfg");

	std::cout<<"Creating rhs"<<std::endl;
	SGPP::base::DataVector b(gridsize);
	b.setAll(0.0);
	operation_mult->generateb(dataset,b);
	for(size_t i=0;i<300;i++)
		std::cout<<b[i]<<" ";
	std::cout<<std::endl;

	std::cout<<"Creating alpha"<<std::endl;
	alpha.setAll(1.0);
	for(unsigned int i=0;i<gridsize;i++)
		alpha[i]=(rand()%1000)/1000.0;
	double norm=alpha.l2Norm();
	alpha.mult(1.0/norm);
	solver->solve(*operation_mult,alpha,b,false,true);

	std::cout<<"Starting graph creation..."<<std::endl;
	SGPP::datadriven::StreamingOCLMultiPlatform::OperationCreateGraphOCL* operation_graph=
		SGPP::datadriven::createNearestNeighborGraphConfigured(dataset, k, dimension, "MyOCLConf.cfg");
	std::vector<int> graph(dataset.getNrows()*k);
	operation_graph->create_graph(graph);

	std::cout<<"Starting graph pruning"<<std::endl;
	SGPP::datadriven::StreamingOCLMultiPlatform::OperationPruneGraphOCL* operation_prune=
		SGPP::datadriven::pruneNearestNeighborGraphConfigured(*grid, dimension, alpha, dataset, treshold, k, "MyOCLConf.cfg");
	operation_prune->prune_graph(graph);

	SGPP::datadriven::StreamingOCLMultiPlatform::OperationCreateGraphOCL::find_clusters(graph, k);

	//cleanup
	delete operation_mult;
	delete solver;


}
