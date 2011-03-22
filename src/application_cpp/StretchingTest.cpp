/*
#include "sgpp.hpp"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstddef>

int main (){
	//write tests for the f-ing strecher :)
	size_t dim =2;
	char c;
	double strike =1;
	sg::Stretching1D* stretching1Ds = new sg::Stretching1D[dim];
	sg::DimensionBoundary* BoundaryArray = new sg::DimensionBoundary[dim];

	for(size_t i =0;i<dim;i++){
		BoundaryArray[i].leftBoundary = 0.1;
		BoundaryArray[i].rightBoundary = 5;
		BoundaryArray[i].bDirichletLeft = true;
		BoundaryArray[i].bDirichletRight = true;
	}
	//	std::string* type = new std::string[dim];
	std::string s0("none");
	std::string s1("log");
	std::string s2("sinh");

	for(size_t j=0;j<dim;j++){
		stretching1Ds[j].type.assign(s0);
	}

	//	stretching1Ds[1].type.assign(s2);
	//	stretching1Ds[1].x_0=3.0;
	//	stretching1Ds[1].xsi=0.75;

	//	std::cout<<type[0]<<std::endl<<type[1]<<std<::endl;
	//	std::cout<<s0<<std::endl<<s1<<std<::endl;

	sg::Stretching strech(dim,BoundaryArray,stretching1Ds);
	sg::Grid* myGrid;
	myGrid = new sg::LinearStretchedTrapezoidBoundaryGrid(strech);

	sg::GridGenerator* myGenerator = myGrid->createGridGenerator();
	myGenerator->regular(4);
	delete myGenerator;

	sg::GridStorage* myGridStorage;
	sg::Stretching* myStretching;
	myStretching = myGrid->getStretching();

	myGridStorage = myGrid->getStorage();
	std::cout<<"Storage Size: "<<myGridStorage->size()<<std::endl;
	DataVector* alpha = new DataVector(myGridStorage->size());
	double tmp;

	//	sg::BoundingBox BB(dim,BoundaryArray);
	//	sg::GridStorage* myGridStorage2;
	//	myGridStorage2 = new sg::GridStorage(BB);
	////	std::string coords2 = myGridStorage2->get(1)->getCoordsStringBB(*(myGridStorage2->getBoundingBox()));
	////	std::cout<<coords2;
	std::cout<<"GetCoordinates: "<<strech.getCoordinates(1,3,1)<<std::endl;


	sg::BoundingBox BB(dim,BoundaryArray);
		sg::Grid* myGrid2;
		myGrid2 = new sg::LinearTrapezoidBoundaryGrid(BB);

		sg::GridGenerator* myGenerator2 = myGrid2->createGridGenerator();
		myGenerator2->regular(11);
		delete myGenerator2;

		sg::GridStorage* myGridStorage2;
		sg::BoundingBox* myBB;
		myBB = myGrid2->getBoundingBox();

		myGridStorage2 = myGrid->getStorage();
		std::cout<<"Storage Size (BB) : "<<myGridStorage2->size()<<std::endl;
		DataVector* alpha2 = new DataVector(myGridStorage2->size());


	for (size_t i = 0; i < dim; i++)
	{
		std::string coords = myGridStorage->get(i)->getCoordsStringStretching(*(myGrid->getStretching()));
		std::cout<<"GetCoordsStretching: "<<coords<<endl;
		std::stringstream coordsStream(coords);
		double* dblFuncValues = new double[dim];

		for (size_t j = 0; j < dim; j++)
		{
			coordsStream >> tmp;
			//std::cout<<tmp<<endl;
			dblFuncValues[j] = tmp;
		}

		tmp = 0.0;
		for (size_t j = 0; j < dim; j++)
		{
			tmp += dblFuncValues[j];
		}

//		std::cout<<"i: "<<i<<endl;

		//std::cout<<alpha[i]<<endl;
		double alphatemp = std::max<double>(((tmp/static_cast<double>(dim))-strike), 0.0);
		alpha->set(i,alphatemp);
//		std::cout<<"alphatemp: "<<alphatemp<<endl;

//		alpha[i] = std::max<double>(((tmp/static_cast<double>(dim))-strike), 0.0);

		std::cout<<"Second tmp Print: "<<tmp<<endl;


		delete[] dblFuncValues;
	}


	for (size_t i = 0; i < dim; i++)
		{
			std::string coords2 = myGridStorage2->get(i)->getCoordsStringBB(*(myGrid2->getBoundingBox()));
			std::cout<<"GetCoordsBB: "<<coords2<<endl;
			std::stringstream coordsStream2(coords2);
			double* dblFuncValues2 = new double[dim];

			for (size_t j = 0; j < dim; j++)
			{
				coordsStream2 >> tmp;
				dblFuncValues2[j] = tmp;
			}

			tmp = 0.0;
			for (size_t j = 0; j < dim; j++)
			{
				tmp += dblFuncValues2[j];
			}

			double alphatemp2 = std::max<double>(((tmp/static_cast<double>(dim))-strike), 0.0);
			alpha2->set(i,alphatemp2);

			delete[] dblFuncValues2;
		}
		sg::OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
		myHierarchisation->doHierarchisation(*alpha);
		delete myHierarchisation;
		std::vector<double> point = vector<double>(2,0.3);
		std::cout<<"Vector values : "<<point[0]<<" "<<point[1]<<endl;

		sg::OperationEval* myEvaluation = myGrid->createOperationEval();
//		sg::OperationEvalLinearStretchedBoundary* myEvaluation =
//				dynamic_cast<sg::OperationEvalLinearStretchedBoundary*>(myGrid->createOperationEval());
		std::cout<<"evaluate (0.3,0.3): "<<myEvaluation->eval(*alpha, point)<<endl;
		delete myEvaluation;


		sg::OperationHierarchisation* myHierarchisation2 = myGrid2->createOperationHierarchisation();
			myHierarchisation2->doHierarchisation(*alpha2);
			delete myHierarchisation2;

			sg::OperationEval* myEvaluation2 = myGrid2->createOperationEval();
//			sg::OperationEvalLinearBoundary* myEvaluation2 =
//					dynamic_cast<sg::OperationEvalLinearBoundary*>(myGrid2->createOperationEval());
			std::cout<<"evaluate (0.3,0.3) for BB : "<<myEvaluation2->eval(*alpha, point)<<endl;
			delete myEvaluation2;

//	std::cout<<strech.getCoordinates(1,3,1)<<std::endl;
	std::cout<<"program is about to end, press a button\n";
	//	std::cin >> &c;

//	delete myGridStorage;
//	delete myStretching;
	delete myGrid;
	delete[] BoundaryArray;
	delete[] stretching1Ds;
	delete alpha;


	//	return 0;
}
 */

/*
#include <iostream>
#include <vector>
// All SG++ headers
#include "sgpp.hpp"
// Or, include only those that are required
//#include "data/DataVector.hpp"
//#include "grid/Grid.hpp"
//#include "grid/GridStorage.hpp"
//#include "grid/generation/GridGenerator.hpp"
//#include "operation/common/OperationEval.hpp"

using namespace std;
using namespace sg;

// function to interpolate
double f(double x0, double x1) {
	return 16.0 * (x0-1)*x0 * (x1-1)*x1;
}

int main() {
	// create a two-dimensional piecewise bi-linear grid
	size_t dim = 2;
	int level =5;
	Grid* grid = Grid::createLinearTrapezoidBoundaryGrid(dim);
	GridStorage* gridStorage = grid->getStorage();
	//	cout << "dimensionality:         " << gridStorage->dim() << endl;

	// create regular grid, level 3
	GridGenerator* gridGen = grid->createGridGenerator();
	gridGen->regular(level);
	//	cout << "number of grid points:  " << gridStorage->size() << endl;

	// create coefficient vector
	DataVector alpha(gridStorage->size());
	alpha.setAll(0.0);
	//	cout << "length of alpha-vector: " << alpha.getSize() << endl;

	// set function values in alpha
	GridIndex* gp;
	for (size_t i=0; i < gridStorage->size(); i++) {
		gp = gridStorage->get(i);
		alpha[i] = f(gp->abs(0), gp->abs(1));
	}
	//	cout << alpha.toString() << endl;

	// hierarchize
	grid->createOperationHierarchisation()->doHierarchisation(alpha);
		cout << alpha.toString() << endl;

	// evaluate
	DataVector p(dim);
	p[0] = 0.52;
	p[1] = 0.73;
	OperationEval* opEval = grid->createOperationEval();
		cout << "For the trivial grid\nu(0.52, 0.73) = " << opEval->eval(alpha, p) << endl;


	DimensionBoundary* BoundaryArray = new DimensionBoundary[dim];

	for(size_t i =0;i<dim;i++){
		BoundaryArray[i].leftBoundary = 0.1;
		BoundaryArray[i].rightBoundary = 2;
		BoundaryArray[i].bDirichletLeft = true;
		BoundaryArray[i].bDirichletRight = true;
	}
	Stretching1D* stretching1Ds = new Stretching1D[dim];
	string s0("none");
	string s1("log");
	string s2("sinh");

	for(size_t j=0;j<dim;j++){
		stretching1Ds[j].type.assign(s2);
		stretching1Ds[j].x_0 = 0.5;
		stretching1Ds[j].xsi=1.0;
	}

	Stretching* strech = new Stretching(dim,BoundaryArray,stretching1Ds);

	Grid* stretchedGrid = new LinearStretchedTrapezoidBoundaryGrid(*strech);
	GridStorage* stretchedGridStorage = stretchedGrid->getStorage();
	//	cout << "stretchedgridStorage dimensionality:         " << stretchedGridStorage->dim() << endl;

	// create regular grid, level 3
	GridGenerator* stretchedGridGen = stretchedGrid->createGridGenerator();
	stretchedGridGen->regular(level);
	//	cout << "number of grid points:  " << stretchedGridStorage->size() << endl;

	// create coefficient vector
	DataVector stretchedAlpha(stretchedGridStorage->size());
	stretchedAlpha.setAll(0.0);
	//	cout << "length of alpha-vector: " << stretchedAlpha.getSize() << endl;

	// set function values in alpha
	GridIndex* sgp;
	for (size_t i=0; i < stretchedGridStorage->size(); i++) {
		sgp = stretchedGridStorage->get(i);
		stretchedAlpha[i] = f(sgp->getCoordStretching(0,strech), sgp->getCoordStretching(1,strech));
	}
	//	cout <<"stretchedAlpha vector: " << stretchedAlpha.toString() << endl;

	// hierarchize
	stretchedGrid->createOperationHierarchisation()->doHierarchisation(stretchedAlpha);
//	cout <<"stAlpha vector after Hierarchization" << stretchedAlpha.toString() << endl;

	// evaluate
	DataVector p2(dim);
	p2[0] = 0.52;
	p2[1] = 0.73;
	OperationEval* opEval2 = stretchedGrid->createOperationEval();
	cout << "For the Stretched Grid\nu(0.52, 0.73) = " << opEval2->eval(stretchedAlpha, p2) << endl;

	double hl=0, hr=0, currentPosition=0;
	//	strech->getAdjacentPositions(4,3,0,hl,hr );
	////currentPosition = strech->getCoordinates(static_cast<int>(current_level), static_cast<int>(current_index), dim);
	//	currentPosition = strech->getCoordinates(4, 3, 0);
	//	cout<<"hl:"<<hl<<" hr:"<<hr<<" curr pos:"<<currentPosition<<endl;
	//strech->printLookupTable();

	BoundingBox BB(dim,BoundaryArray);

	Grid* BBGrid = new LinearTrapezoidBoundaryGrid(BB);
	GridStorage* BBGridStorage = BBGrid->getStorage();
	//	cout << "stretchedgridStorage dimensionality:         " << BBGridStorage->dim() << endl;

	// create regular grid, level 3
	GridGenerator* BBGridGen = BBGrid->createGridGenerator();
	BBGridGen->regular(level);
	//	cout << "number of grid points:  " << BBGridStorage->size() << endl;

	// create coefficient vector
	DataVector BBAlpha(BBGridStorage->size());
	BBAlpha.setAll(0.0);
	//	cout << "length of alpha-vector: " << BBAlpha.getSize() << endl;

	// set function values in alpha
	GridIndex* bbgp;
	double q0 = BB.getIntervalWidth(0);
	double t0 = BB.getIntervalOffset(0);
	double q1 = BB.getIntervalWidth(1);
		double t1 = BB.getIntervalOffset(1);
	for (size_t i=0; i < BBGridStorage->size(); i++) {
		bbgp = BBGridStorage->get(i);
		BBAlpha[i] = f(bbgp->getCoordBB(0,q0,t0), bbgp->getCoordBB(1,q1,t1));
	}
	//	cout <<"BBAlpha vector: " <<BBAlpha.toString() << endl;

	// hierarchize
	BBGrid->createOperationHierarchisation()->doHierarchisation(BBAlpha);
//	cout << "BBAlpha vector after Hierarchization"<<BBAlpha.toString() << endl;

	// evaluate
	DataVector p3(dim);
	p3[0] = 0.52;
	p3[1] = 0.73;
	OperationEval* opEval3 = BBGrid->createOperationEval();
	cout << "For the BB Grid\nu(0.52, 0.73) = " << opEval3->eval(BBAlpha, p3) << endl;



 * Part where the testing of Janos Stretching is done.

	double myints[] = {0, 0.5, 1};
	vector<double> v1(myints, myints + sizeof(myints) / sizeof(double) );
	vector<double>* myVecArr = new vector<double>[dim];
	myVecArr[0]=v1;
	myVecArr[1]=v1;

	cout<<myVecArr[0][0]<<myVecArr[1][4]<<endl;
	Stretching strech2(dim, myVecArr);
	strech2.printLookupTable();




//	delete [] BoundaryArray;
	delete [] stretching1Ds;
//	delete BBGrid;
////	delete grid;
	delete strech;
	delete stretchedGrid;


	return 0;

}
 */


#include <iostream>
#include <vector>
#include "ctime"
// All SG++ headers
#include "sgpp.hpp"

using namespace std;
using namespace sg;

// function to interpolate
double f(double x0, double x1) {
	return 16.0 * (x0-1)*x0 * (x1-1)*x1;
}

void waitt ( int seconds )
{
  clock_t endwait;
  endwait = clock () + seconds * CLOCKS_PER_SEC ;
  while (clock() < endwait) {}
}
int main() {
	// create a two-dimensional piecewise bi-linear grid
	size_t dim = 2;
	int level =3;

	DimensionBoundary* BoundaryArray = new DimensionBoundary[dim];

	for(size_t i =0;i<dim;i++){
		BoundaryArray[i].leftBoundary = 0.5;
		BoundaryArray[i].rightBoundary = 7;
		BoundaryArray[i].bDirichletLeft = true;
		BoundaryArray[i].bDirichletRight = true;
	}
	Stretching1D* stretching1Ds = new Stretching1D[dim];
	string s0("none");
	string s1("log");
	string s2("sinh");

	for(size_t j=0;j<dim;j++){
		stretching1Ds[j].type.assign(s1);
		stretching1Ds[j].x_0 = 1.0;
		stretching1Ds[j].xsi=10;
	}
	//	stretching1Ds[0].type.assign(s0);

	Stretching* strech = new Stretching(dim,BoundaryArray,stretching1Ds);

	//	strech->printLookupTable();



	Grid* stretchedGrid = new LinearStretchedTrapezoidBoundaryGrid(dim);
	stretchedGrid->setStretching(*strech);

	//	cout << "stretchedgridStorage dimensionality:         " << stretchedGridStorage->dim() << endl;

	// create regular grid, level 3
	GridGenerator* stretchedGridGen = stretchedGrid->createGridGenerator();
	stretchedGridGen->regular(level);
	//	cout << "number of grid points:  " << stretchedGridStorage->size() << endl;
	GridStorage* stretchedGridStorage = stretchedGrid->getStorage();

	// create coefficient vector
	DataVector stretchedAlpha(stretchedGrid->getStorage()->size());
	//	DataVector stretchedAlpha(stretchedGridStorage->size());
	stretchedAlpha.setAll(0.0);
	//		cout << "length of alpha-vector: " << stretchedAlpha.getSize() << endl;

	////	 set function values in alpha
	//	GridIndex* sgp;
	//	for (size_t i=0; i < stretchedGrid->getStorage()->size(); i++) {
	//
	//		sgp = stretchedGrid->getStorage()->get(i);
	//		stretchedAlpha[i] = f(sgp->getCoordStretching(0,stretchedGrid->getStorage()->getStretching()), sgp->getCoordStretching(1,stretchedGrid->getStorage()->getStretching()));
	////		stretchedAlpha[i] = f(sgp->getCoordStretching(0,stretchedGrid->getStorage()->getStretching()), sgp->getCoordStretching(0,stretchedGrid->getStorage()->getStretching()));

	//}
	//	cout <<"stretchedAlpha vector: " << stretchedAlpha.toString() << endl;

	// hierarchize
	stretchedGrid->createOperationHierarchisation()->doHierarchisation(stretchedAlpha);
	//	cout <<"stAlpha vector after Hierarchization" << stretchedAlpha.toString() << endl;
	//
	//	// evaluate
	//	DataVector p2(dim);
	//	p2[0] = 0.52;
	//	p2[1] = 0.73;
	//	OperationEval* opEval2 = stretchedGrid->createOperationEval();
	//	cout << "For the Stretched Grid\nu(0.52, 0.73) = " << opEval2->eval(stretchedAlpha, p2) << endl;
	//
	//	double hl=0, hr=0, currentPosition=0;

	/*	string storageString = stretchedGridStorage->serialize();
	cout<<storageString<<endl;
//
//
//
//	BoundingBox BB(dim,BoundaryArray);
//
//		Grid* BBGrid = new LinearTrapezoidBoundaryGrid(BB);
//		GridStorage* BBGridStorage = BBGrid->getStorage();
//
//		string storageString2 = BBGridStorage->serialize();
//			cout<<storageString2<<endl;
//
//			GridStorage* newstretchedGridStorage = new GridStorage(storageString);
//			string storageString3 = newstretchedGridStorage->serialize();
//				cout<<storageString3<<endl;*/
	////	std::cout<<"Print Sparse Grid\n";
	////	GridPrinterForStretching gpfs= GridPrinterForStretching(*stretchedGrid);
	////
	////		gpfs.printSparseGrid(stretchedAlpha, "gridi_print.txt", false);

	//	cout<< "Test getDiscreteVector\n";
	double myDoubles[] = {0,1,2,3,4,5,6,7,8};
	vector<double> vec1(myDoubles,  myDoubles + sizeof(myDoubles) / sizeof(double) );

	vector<double>* vec = new vector<double>[dim];
	vector<double>* vecc;
	for(size_t i=0;i<dim;i++){
		vec[i]=vec1;
	}

	Stretching* strechDisc = new Stretching(dim, vec);

	double myDoubles2[] =  {1.05,0.553045,1.58554,0.29184,0.80243,1.30366,2.03825,0.168726,	0.423203,0.678957,0.925664,	1.1758,	1.43689,
			1.7721,
			2.43901};
	vector<double> vec2(myDoubles2,  myDoubles2 + sizeof(myDoubles2) / sizeof(double) );

	string* osman = strechDisc->getStretchingMode();

//	for(int i=0; i<15;i++){
//	cout<<vec2[i]<<" ";
//	}
//	cout<<endl;

	vecc = strechDisc -> getDiscreteVector(true);
	//for(size_t i=0; i<dim; i++){
	//	for(int j=0; j<vecc[i].size();j++ ){
	//		 cout << " " << vecc[i][j];
	//	}
	//	cout<<endl;
	//}

/*
 * Print the serialization
 */
/*	Grid* stretchedGrid2 = new LinearStretchedTrapezoidBoundaryGrid(dim);
	stretchedGrid2->setStretching(*strechDisc);
	GridGenerator* stretchedGridGen2 = stretchedGrid2->createGridGenerator();
	stretchedGridGen2->regular(level);
	string ebuuve(stretchedGrid2->serialize());
	cout<<ebuuve;

	GridPrinterForStretching gpfs2= GridPrinterForStretching(*stretchedGrid2);

	gpfs2.printSparseGrid(stretchedAlpha, "grid2.txt", false);

	Grid* stretchedGrid3 = stretchedGrid2->unserialize(ebuuve);
	cout<<stretchedGrid3->serialize();

	GridPrinterForStretching gpfs3= GridPrinterForStretching(*stretchedGrid3);
	gpfs3.printSparseGrid(stretchedAlpha, "grid3.txt", false);*/

//	cout<<"Get Coordinates:\n"<<strechDisc->getCoordinates(15,501,1)<<endl;

//	Grid* stretchedGrid2 = new LinearStretchedTrapezoidBoundaryGrid(dim);
//		stretchedGrid2->setStretching(*strechDisc);
//		GridGenerator* stretchedGridGen2 = stretchedGrid2->createGridGenerator();
//		stretchedGridGen2->regular(level);
//		string ebuuve(stretchedGrid2->serialize());
//		cout<<ebuuve;

//		Grid* SGpp_grid_ = sg::Grid::createLinearStretchedTrapezoidBoundaryGrid(dim);
//		// set stretching
//		sg::Stretching* myStretching = strechDisc;
//		SGpp_grid_->setStretching(*myStretching);
//
////		GridGenerator* stretchedGridGen3 = SGpp_grid_->createGridGenerator();
////		SGpp_grid_->regular(level);
//		std::string serializedGrid = "";
//		SGpp_grid_->serialize(serializedGrid);
//		cout <<serializedGrid;
//
//
//		BlackScholesSolver* myBSSolver_;
//bool withStretching = true;
//		if(withStretching){
//		 myBSSolver_ = new BlackScholesSolverWithStretching(false, "European");
//		}else{
//		  myBSSolver_ = new BlackScholesSolver( false , "European");
//		}
//		 myBSSolver_->initScreen();
//		 myBSSolver_->constructGrid(*myStretching,3);
//
//		 cout << typeid(*myBSSolver_).name() << endl;
//
//		 BlackScholesSolverWithStretching* myBSSolver2_;
//		 myBSSolver2_ = new sg::BlackScholesSolverWithStretching(
//		     false, "European");
//		 myBSSolver2_->initScreen();
	double posr,posl,posc;
	int maxindex=0;

time_t istart=0,iend=0;

clock_t cstart=clock();
time(&istart);
for(int counter=0; counter<1000;counter++){
	for(size_t d=0;d<dim; d++){
		for(int l=0;l<15;l++){
			maxindex=pow(2,l);
			for(int i=1;i<maxindex;i=i+2){
//				cout<<"\n";
				strech->getAdjacentPositions(l,i,d,posc, posl,posr);
//				posc=i;
			}
		}
	}
}
	time(&iend);
	clock_t  cend = clock();
cout<<scientific<<"istart "<<istart<<endl;
cout<<scientific<<"end "<<iend<<endl;

cout<<scientific<<"cstart "<<cstart<<endl;
cout<<scientific<<"cend "<<cend<<endl;

cout<<"Takes "<<difftime(iend, istart)<<" seconds\n";
cout<<"Takes "<<double(cend-cstart)/CLOCKS_PER_SEC<<" seconds\n";

clock_t t1=clock();
for(int i=0; i<1000000;i++);
clock_t t2=clock();
cout<<"Iterations took for " << double(t2-t1)/CLOCKS_PER_SEC << " seconds\n";


	delete [] BoundaryArray;
	delete [] stretching1Ds;
	delete strech;
	delete stretchedGrid;


//  int n;
//
//  clock_t start;
//
//  clock_t end;
//
//  int duration;
//
//  start = clock();
//
//  printf ("Starting countdown...\n");
//
//  for (n=10; n>0; n--)
//
//  {
//
//    printf ("%d\n",n);
//
//    waitt (1);
//
//  }
//
//  end = clock();
//
//  duration = (int) (end-start);
//
//  printf ("FIRE!!! %d\n", duration);


	return 0;

}



