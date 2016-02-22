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
	size_t dimension,tiefe,k;
	double lambda,treshold;
	std::string filename = "dataset2_dim2.arff";

	std::cout << "Loading file: " << filename << std::endl;
	SGPP::datadriven::Dataset data =
		SGPP::datadriven::ARFFTools::readARFF(filename);
	SGPP::base::DataMatrix& dataset = data.getData();
	dimension=dataset.getNcols();
	std::cout<<"Loaded "<<dataset.getNcols()<<" dimensional dataset with "<<dataset.getNrows()<<" datapoints."<<std::endl;

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

	sg::solver::ConjugateGradients *solver=new sg::solver::ConjugateGradients(1000,0.001);
	SGPP::datadriven::StreamingOCLMultiPlatform::OperationDensityOCL* operation_mult=
		SGPP::datadriven::createDensityOCLMultiPlatformConfigured(*grid, dimension, lambda, "MyOCLConf.cfg");

	std::cout<<"Creating rhs"<<std::endl;
	SGPP::base::DataVector b(gridsize);
	operation_mult->generateb(dataset,b);
	for(size_t i=0;i<300;i++)
		std::cout<<b[i]<<" ";
	std::cout<<std::endl;

	std::cout<<"Creating alpha"<<std::endl;
	solver->solve(*operation_mult,alpha,b,false,true);
	double max=alpha.max();
	double min=alpha.min();
	for(size_t i=0;i<gridsize;i++)
		alpha[i]=alpha[i]*1.0/(max-min);

	std::cout<<"Starting graph creation..."<<std::endl;
	SGPP::datadriven::StreamingOCLMultiPlatform::OperationCreateGraphOCL* operation_graph=
		SGPP::datadriven::createNearestNeighborGraphConfigured(dataset, k, dimension, "MyOCLConf.cfg");
	std::vector<int> graph(dataset.getNrows()*k);
	operation_graph->create_graph(graph);
	std::ofstream out1("graph1.txt");
	for(size_t i=0; i < dataset.getNrows(); i++) {
		for(size_t node = 0; node < k; node++)
			out1 << graph[i * k + node] << " ";
		out1 << std::endl;
	}
	out1.close();

	std::cout<<"Starting graph pruning"<<std::endl;
	SGPP::datadriven::StreamingOCLMultiPlatform::OperationPruneGraphOCL* operation_prune=
		SGPP::datadriven::pruneNearestNeighborGraphConfigured(*grid, dimension, alpha, dataset, treshold, k, "MyOCLConf.cfg");
	operation_prune->prune_graph(graph);

	std::ofstream out2("graph2.txt");
	for(size_t i=0; i < dataset.getNrows(); i++) {
		for(size_t node = 0; node < k; node++)
			out2 << graph[i * k + node] << " ";
		out2 << std::endl;
	}
	out2.close();

	SGPP::datadriven::StreamingOCLMultiPlatform::OperationCreateGraphOCL::find_clusters(graph, k);

	//cleanup
	delete gridGen;
	delete grid;
	delete operation_mult;
	delete solver;


}
