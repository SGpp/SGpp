#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/exception/file_exception.hpp>
#include <sgpp/base/exception/generation_exception.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/exception/solver_exception.hpp>
#include <sgpp/base/exception/tool_exception.hpp>

//For algorithm exception
#include <sgpp/base/datatypes/DataVector.hpp>

//For application exception
//#include <sgpp/datadriven/datamining/DataMiningConfiguration.hpp>
//#include <string>

//For data exception
#include <sgpp/base/datatypes/DataMatrix.hpp>

//For factory exception
#include <sgpp/base/tools/QuadRule1D.hpp>

//For file exception
#include <sgpp/base/grid/GridDataBase.hpp>

//For generic exception
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/generation/PeriodicGridGenerator.hpp>

//For operation exception
#include <sgpp/base/operation/BaseOpFactory.hpp>

//For solver exception
//#include <sgpp/solver/ode/Euler.hpp>
//#include <sgpp/base/application/ScreenOutput.hpp>

//For tool exception
#include <sgpp/base/tools/GridPrinter.hpp>

#include <sgpp/globaldef.hpp>

using namespace SGPP::base;
//using namespace SGPP::solver;


BOOST_AUTO_TEST_CASE(test_AlgorithmException){
	DataVector testVector(5);
	testVector.setAll(0);
	algorithm_exception expectedException("more indices than entries!");
	std::vector<size_t> tooLongVector(6);
	try{
		testVector.restructure(tooLongVector);
	}catch(algorithm_exception& actualException){
		BOOST_CHECK_EQUAL(actualException.what(), expectedException.what());
	}
}

//TODO: Add this when datadriven is working
//BOOST_AUTO_TEST_CASE(test_ApplicationException){
//	SGPP::datadriven::DataMiningConfiguration* dataMiningConfigurationDummy;
//	std::string str = "haha... this is not a grid type!";
//	application_exception expectedException("grid type is unknown");
//
//	try{
//		dataMiningConfigurationDummy->stringToGridType(str);
//	}catch(application_exception& actualException){
//		BOOST_CHECK_EQUAL(actualException.what(), expectedException.what());
//	}
//}

BOOST_AUTO_TEST_CASE(test_DataException){
	data_exception expectedException("DataVector::add : Dimensions do not match");
	DataVector dataVectorOne(3, 1.0);
	DataVector dataVectorTwo(4, 2.9);

	try{
		dataVectorOne.add(dataVectorTwo);
	}catch(data_exception* actualException){
		BOOST_CHECK_EQUAL(actualException->what(), expectedException.what());
	}
}

BOOST_AUTO_TEST_CASE(test_FactoryException){
	QuadRule1D quadRule1D;

	size_t faultyLevel = -1;
	DataVector pweight(5);
	DataVector coordinates(5);
    factory_exception expectedException(
      "QuadRule1D::getLevelPointsAndWeights : "
      "order of gauss quadrature has to be within {1, ..., 20}");
	try{
		quadRule1D.getLevelPointsAndWeights(faultyLevel, coordinates, pweight);
	}catch(factory_exception& actualException){
		BOOST_CHECK_EQUAL(actualException.what(), expectedException.what());
	}
}

BOOST_AUTO_TEST_CASE(test_FileException){
	size_t dimension = 2;
	GridDataBase gridDataBase(dimension);
	std::string filename = "fooFooFOOFile1234";
	std::string expectedMessage = "Error! Unable to open file '" + filename + "' for read access.";
	file_exception expectedException(expectedMessage.c_str());
	try{
		gridDataBase.save(filename);
	}catch(file_exception* actualException){
		BOOST_CHECK_EQUAL(actualException->what(), expectedException.what());
	}
	std::remove(filename.c_str());
}

BOOST_AUTO_TEST_CASE(test_GenerationException){
	size_t dimension = 2;
	Grid* linearGrid = Grid::createLinearGrid(dimension);
	GridStorage* linearGridStorage = linearGrid->getStorage();

	PeriodicGridGenerator gridGen(linearGridStorage);
	size_t level = 3;
	generation_exception expectedException("PeriodicGridGenerator::full is not implemented");

	try{
		gridGen.full(level);
	}catch(generation_exception& actualException){
		BOOST_CHECK_EQUAL(actualException.what(), expectedException.what());
	}
}

BOOST_AUTO_TEST_CASE(test_OperationException){
	size_t dimension = 7;
	Grid* linearGrid = Grid::createLinearGrid(dimension);
	DataMatrix dummyMatrix(10,7);
	OperationMultipleEval* opEval = SGPP::op_factory::createOperationMultipleEval(*linearGrid, dummyMatrix);

	operation_exception expectedException(
		      "error: OperationMultipleEval::getDuration(): "
		      "not implemented for this kernel");

	try{
		opEval->getDuration();
	}catch(operation_exception* actualException){
		BOOST_CHECK_EQUAL(actualException->what(), expectedException.what());
	}
}

//TODO: Got a linking error here. Not sure why!
//BOOST_AUTO_TEST_CASE(test_SolverException){
//    std::string faultyString = "haha... I am not an Euler mode";
//	size_t imax = 10;
//	float_t timestepsize = 0.1f;
//	bool generateAnimation = false;
//	size_t numEvalsAnimation = 10;
//	ScreenOutput screen;
//	solver_exception expectedException("Euler::Euler : An unknown Euler-Mode was specified!");
//	Euler* faultyEuler;
//	try{
//		faultyEuler = new SGPP::solver::Euler(faultyString, imax, timestepsize, generateAnimation, numEvalsAnimation, &screen);
//	}catch(solver_exception* actualException){
//		BOOST_CHECK_EQUAL(actualException->what(), expectedException.what());
//	}
//	delete faultyEuler;
//}

BOOST_AUTO_TEST_CASE(test_ToolException){
	size_t dimension = 4;
	Grid* linearGrid = Grid::createLinearGrid(dimension);
	GridPrinter gridPrinter(*linearGrid);

	tool_exception expectedException("GridPrinter::printGrid : "
            "The grid has no dimensions. "
            "Thus it cannot be printed!");
	size_t pointsPerDim = 10;
	std::string dummyFilename = "fooFooFooef";
	DataVector dummyVector(10);
	try{
		gridPrinter.printGrid(dummyVector, dummyFilename, pointsPerDim);
	}catch(tool_exception* actualException){
		BOOST_CHECK_EQUAL(actualException->what(), expectedException.what());
	}
}
