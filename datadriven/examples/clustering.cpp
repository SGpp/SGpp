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

using namespace SGPP::base;
int main()
{
	unsigned int dimension,tiefe,k;
	double lambda,treshold;
	std::string filename;
	std::cout<<"Count of dimensions: ";
	std::cin>>dimension;
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
	//Load dataset
	filename="simple_test.txt";
	std::ifstream teststream(filename.c_str(),std::ifstream::in);
	char zeile[256];
	teststream.getline(zeile,256);
	std::string zeilenstring(zeile);
	std::stringstream stream(zeilenstring);
	int spaltencounter=0;
	int zeilencounter=1;
	double dummy=0;
	while(stream.good())
	{
		stream>>dummy;
		spaltencounter++;
	}
	while(teststream.getline(zeile,256))
		zeilencounter++;
	std::cout<<"Number of datapoints:"<<zeilencounter<<std::endl;
	teststream.close();
	std::ifstream ifs(filename.c_str(),std::ifstream::in);
	float tmp=0;
	unsigned int counter=0;
	std::vector<float> buffer_data;
	while(counter<zeilencounter*spaltencounter)
	{
		ifs>>tmp;
		buffer_data.push_back(tmp);
		counter++;
	}
	ifs.close();
	std::cout<<"Dataset loaded!"<<std::endl;

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
	SGPP::base::DataVector dataset(buffer_data.size());
	for(size_t i=0; i < buffer_data.size();i++)
		dataset[i]=buffer_data[i];
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
	std::vector<int> graph(zeilencounter*k);
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
