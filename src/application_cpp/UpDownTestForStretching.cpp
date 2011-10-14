/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include "sgpp.hpp"
#include "exception/operation_exception.hpp"
#include "grid/type/LinearStretchedGrid.hpp"
//#include "grid/type/LinearBoundaryGrid.hpp"
#include "grid/type/LinearStretchedTrapezoidBoundaryGrid.hpp"
//#include "grid/type/ModLinearGrid.hpp"

#include <iostream>
#include <string>
#include <fstream>

#define PRINTMATRIX
//#define PRINTDIFF
//#define BOUNDING

int main(int argc, char *argv[])
{
	size_t levels = 5;
	size_t dim = 3;
	size_t numGridPoints;

	sg::base::Grid* myGrid1;
	sg::base::OperationMatrix* myUpDown1;
	sg::base::OperationMatrix* myDelta1;
	sg::base::OperationMatrix* myGamma1;
	sg::base::OperationMatrix* myLaplace1;
	sg::base::OperationMatrix* myL2Dot1;

	std::cout << std::endl;
	std::cout << "Starting the Up / Down Test" << std::endl;
	std::cout << "===========================" << std::endl;
	std::cout << std::endl;
	std::cout << "levels:" << levels << std::endl;
	std::cout << "dim:   " << dim << std::endl;
	std::cout << std::endl;


	sg::base::Stretching* myStretching;
	sg::base::DimensionBoundary* myInterval;

	if(dim==1){

		myInterval = new sg::base::DimensionBoundary;
		myInterval->leftBoundary = 0.5;
		myInterval->rightBoundary = 7.0;
		myInterval->bDirichletLeft = false;
		myInterval->bDirichletRight = false;
	}
	else if(dim>1)
	{

		myInterval=new sg::base::DimensionBoundary[dim];
		for(size_t j=0; j<dim;j++){
			myInterval[j].leftBoundary = 0.5;
			myInterval[j].rightBoundary = 7.0;
			myInterval[j].bDirichletLeft = false;
			myInterval[j].bDirichletRight = false;
		}
	}


	sg::base::Stretching1D* stretching1Ds = new sg::base::Stretching1D[dim];
	string s0("id");
	string s1("log");
	string s2("sinh");

	for(size_t j=0;j<dim;j++){
		stretching1Ds[j].type.assign(s2);
		stretching1Ds[j].x_0 = 1.0;
		stretching1Ds[j].xsi=10;
	}

	myStretching = new sg::base::Stretching(dim, myInterval, stretching1Ds );
	//	myStretching->printLookupTable();

	std::cout << "start constructing regular grid" << std::endl;
	myGrid1 = new sg::base::LinearStretchedTrapezoidBoundaryGrid((*myStretching));
	std::cout << "A LinearStretchedTrapezoidBoundaryGrid was created" << std::endl;

	sg::base::GridGenerator* myGenerator1 = myGrid1->createGridGenerator();
	myGenerator1->regular(levels);
	delete myGenerator1;
	std::cout << levels << " levels were added to the above created grid" << std::endl;

	numGridPoints = myGrid1->getStorage()->size();
	std::cout << "the grid has " << numGridPoints << " gridpoints" << std::endl;
	std::cout << "finished construction regular grid" << std::endl;
	//std::cout << "the grid is:" << std::endl;

	std::string ser1;
	myGrid1->serialize(ser1);
	//std::cout << ser << std::endl;


	DataVector alpha1(numGridPoints);
	alpha1.setAll(0.0);

	DataVector result1(numGridPoints);
	result1.setAll(0.0);

	DataVector result12(numGridPoints);
	result12.setAll(0.0);

	DataVector result13(numGridPoints);
	result13.setAll(0.0);

	DataVector result14(numGridPoints);
	result14.setAll(0.0);

	DataVector result15(numGridPoints);
	result15.setAll(0.0);

	DataMatrix UpDownMatrix1(numGridPoints, numGridPoints);
	UpDownMatrix1.setAll(0.0);

	DataMatrix GammaMatrix1(numGridPoints, numGridPoints);
	GammaMatrix1.setAll(0.0);

	DataMatrix DeltaMatrix1(numGridPoints, numGridPoints);
	DeltaMatrix1.setAll(0.0);

	DataMatrix LaplaceMatrix1(numGridPoints, numGridPoints);
	LaplaceMatrix1.setAll(0.0);

	DataMatrix L2DotMatrix1(numGridPoints, numGridPoints);
	L2DotMatrix1.setAll(0.0);



	/*	DataVector deltaCoef(numGridPoints),mus(numGridPoints);
	DataMatrix gammaCoef(numGridPoints,numGridPoints);
	gammaCoef.setAll(0.97);
	mus.setAll(0.05);
	DataMatrix rhos(numGridPoints,numGridPoints);
	rhos.setAll(0.0);*/


	DataVector deltaCoef(dim),mus(dim),sigmas(dim);
	DataMatrix gammaCoef(dim,dim);
	gammaCoef.setAll(0);
	DataMatrix rhos(dim,dim);
	rhos.setAll(0.0);


	sigmas.setAll(0.4);
	mus.setAll(0.05);

	for(int k=0; k<dim;k++){
		sigmas[k]=0.4;
		rhos.set(k,k,1.0);
	}

/*	//Update for 2D in Results section
	sigmas[1]=0.5;
	rhos.set(0,1,0.1);
	rhos.set(1,0,0.1);*/

	//Update for 3D in Results section
	sigmas[0] = 0.2; sigmas[1]=0.3; sigmas[2]=0.4;
	mus[0]=0.1; mus[1]=0.02; mus[2]=0.04;
	rhos.set(0,0,1.0); rhos.set(0,1,-0.7); rhos.set(0,2,-0.1);
	rhos.set(1,0,-0.7); rhos.set(1,1,1.0); rhos.set(1,2,0.1);
	rhos.set(2,0,-0.1); rhos.set(2,1,0.1); rhos.set(2,2,1.0);



	for (size_t i = 0; i < dim; i++)
	{
		for (size_t j = 0; j < dim; j++)
		{
			// handle diagonal
			if (i == j)
			{
				gammaCoef.set(i, j, 0.5*(sigmas[i]*sigmas[j]*rhos.get(i,j)));
			}
			else
			{
				gammaCoef.set(i, j, (sigmas[i]*sigmas[j]*rhos.get(i,j)));
			}
		}
	}

	for (size_t i = 0; i < dim; i++)
	{
		for (size_t j = 0; j < dim; j++)
		{
			//					std::cout << sigmas[i] << " ";
			//					std::cout << std::scientific << rhos.get(i, j) << " ";
			//			std::cout << std::scientific << gammaCoef.get(i, j) << " ";
			//					outfile1 << std::scientific << UpDownMatrix1.get(i, j) << " ";
		}
		std::cout << std::endl;
		//				outfile1 << std::endl;
	}

	double covar_sum = 0.0;

	for (size_t i = 0; i < dim; i++)
	{
		covar_sum = 0.0;
		for (size_t j = 0; j < dim; j++)
		{
			// handle diagonal
			if (i == j)
			{
				covar_sum += ((sigmas[i]*sigmas[j]*rhos.get(i,j)));
			}
			else
			{
				covar_sum += (0.5*(sigmas[i]*sigmas[j]*rhos.get(i,j)));
			}
		}
		deltaCoef.set(i, mus[i]-covar_sum);
	}


	myUpDown1 = sg::op_factory::createOperationUpDownTest(*myGrid1);
	myDelta1 = sg::op_factory::createOperationDelta(*myGrid1, deltaCoef);
	myGamma1 = sg::op_factory::createOperationGamma(*myGrid1, gammaCoef);
	myL2Dot1 = sg::op_factory::createOperationLTwoDotProduct(*myGrid1);

	std::cout << "start constructing the operator's matrix" << std::endl;
	for (size_t i = 0; i < numGridPoints; i++)
	{
		// init alpha
		alpha1.setAll(0.0);
		alpha1.set(i, 1.0);
		// init result
		result1.setAll(0.0);
		result12.setAll(0.0);
		result13.setAll(0.0);
//		result14.setAll(0.0);
		result15.setAll(0.0);

		myUpDown1 -> mult(alpha1, result1);
		myDelta1 -> mult(alpha1,result12);
		myGamma1 -> mult(alpha1, result13);
//		myLaplace1-> mult(alpha1, result14);
		myL2Dot1-> mult(alpha1, result15);
		// copy data to operator's matrix
		UpDownMatrix1.setColumn(i, result1);
		DeltaMatrix1.setColumn(i, result12);
		GammaMatrix1.setColumn(i, result13);
//		LaplaceMatrix1.setColumn(i, result14);
		L2DotMatrix1.setColumn(i, result15);
	}
	std::cout << "finished constructing the operator's matrix" << std::endl;

#ifdef PRINTMATRIX
	std::cout << "The operator's matrix is:" << std::endl << std::endl;
	std::ofstream outfile1;


	outfile1.open("updownforstretching.test");

	for (size_t i = 0; i < numGridPoints; i++)
	{
		for (size_t j = 0; j < numGridPoints; j++)
		{
			//			std::cout << std::scientific << UpDownMatrix1.get(i, j) << " ";
			outfile1 << std::scientific << UpDownMatrix1.get(i, j) << " ";
		}
		//		std::cout << std::endl;
		outfile1 << std::endl;
	}

	outfile1.close();

	outfile1.open("deltaforstretching.test");

	for (size_t i = 0; i < numGridPoints; i++)
	{
		for (size_t j = 0; j < numGridPoints; j++)
		{
			//			std::cout << std::scientific << DeltaMatrix1.get(i, j) << " ";
			outfile1 << std::scientific << DeltaMatrix1.get(i, j) << " ";
		}
		//		std::cout << std::endl;
		outfile1 << std::endl;
	}

	outfile1.close();

	outfile1.open("gammaforstretching.test");

	for (size_t i = 0; i < numGridPoints; i++)
	{
		for (size_t j = 0; j < numGridPoints; j++)
		{
			//			std::cout << std::scientific << GammaMatrix1.get(i, j) << " ";
			outfile1 << std::scientific << GammaMatrix1.get(i, j) << " ";
		}
		//		std::cout << std::endl;
		outfile1 << std::endl;
	}

	outfile1.close();


/*	outfile1.open("laplaceforstretching.test");

	for (size_t i = 0; i < numGridPoints; i++)
	{
		for (size_t j = 0; j < numGridPoints; j++)
		{
			//			std::cout << std::scientific << DeltaMatrix1.get(i, j) << " ";
			outfile1 << std::scientific << LaplaceMatrix1.get(i, j) << " ";
		}
		//		std::cout << std::endl;
		outfile1 << std::endl;
	}

	outfile1.close()*/;

	outfile1.open("l2dotproductforstretching.test");

	for (size_t i = 0; i < numGridPoints; i++)
	{
		for (size_t j = 0; j < numGridPoints; j++)
		{
			//			std::cout << std::scientific << DeltaMatrix1.get(i, j) << " ";
			outfile1 << std::scientific << L2DotMatrix1.get(i, j) << " ";
		}
		//		std::cout << std::endl;
		outfile1 << std::endl;
	}

	outfile1.close();
#endif

#ifdef PRINTDIFF
	std::string file1 = "XdPhiPhi1D.txt";
	std::ifstream infile1;

	infile1.open(file1.c_str());

	double filedata1 = 0.0;

	for (size_t i = 0; i < numGridPoints; i++)
	{
		for (size_t j = 0; j < numGridPoints; j++)
		{
			infile1 >> filedata1;
			std::cout << (UpDownMatrix1.get(i, j) - filedata1) << " ";
		}
		std::cout << std::endl;
	}

	infile1.close();
#endif

	std::cout << std::endl;

	std::cout << "Test symmetry:" << std::endl;
	double symTest1 = 0.0;
	double tempSym1 = 0.0;
	for (size_t i = 0; i < numGridPoints; i++)
	{
		for (size_t j = 0; j < numGridPoints; j++)
		{
			tempSym1 = fabs(UpDownMatrix1.get(i,j)-UpDownMatrix1.get(j,i));
			if (tempSym1 > symTest1)
			{
				symTest1 = tempSym1;
			}
		}
	}
	std::cout << "Maximum symmetry error: " << symTest1 << std::endl << std::endl;

	/*	std::cout << "Do a Test multiplication:" << std::endl << std::endl;

	alpha1.setAll(0.0);
	alpha1.set(1, 35.0);
	alpha1.set(2, -17.5);
	alpha1.set(4, -7.5);
	alpha1.set(7, -5);
	alpha1.set(14, -1.25);

	result1.setAll(0.0);

	myUpDown1->mult(alpha1, result1);

	std::cout << result1.toString() << std::endl;

	std::cout << std::endl;*/

	/*sg::IOToolBonnSG* myImporter = new sg::IOToolBonnSG();

	DataVector serAlpha(0);
	std::string serGrid;
	bool hier;

	myImporter->readFile("sparse_bonn.grid", serGrid, serAlpha, hier);

	myImporter->writeFile("sgpp_bonn.grid", *myGrid, alpha, true);

	delete myImporter;*/
	delete myUpDown1;
	delete myGrid1;
#ifdef BOUNDING
	sg::base::Grid* myGrid2;
	sg::base::OperationMatrix* myUpDown2;

	sg::base::BoundingBox* myBoundingBox = new sg::base::BoundingBox(dim, myInterval);

	std::cout << "start constructing regular grid" << std::endl;
	myGrid2 = new sg::LinearTrapezoidBoundaryGrid((*myBoundingBox));
	std::cout << "A LinearTrapezoidBoundaryGrid was created" << std::endl;

	sg::base::GridGenerator* myGenerator2 = myGrid2->createGridGenerator();
	myGenerator2->regular(levels);
	delete myGenerator2;
	std::cout << levels << " levels were added to the above created grid" << std::endl;

	numGridPoints = myGrid2->getStorage()->size();
	std::cout << "the grid has " << numGridPoints << " gridpoints" << std::endl;
	std::cout << "finished construction regular grid" << std::endl;
	//std::cout << "the grid is:" << std::endl;

	std::string ser2;
	myGrid2->serialize(ser2);
	//std::cout << ser << std::endl;



	//		sg::base::DimensionBoundary myInterval;
	//
	//		myInterval.leftBoundary = 0.001;
	//		myInterval.rightBoundary = 100.0;
	//		myInterval.bDirichletLeft = false;
	//		myInterval.bDirichletRight = false;

	//		myBoundingBox->setBoundary(0, myInterval);
	//
	//		std::cout << "Changed Bounding Box to: " << myInterval.leftBoundary << " to " << myInterval.rightBoundary << std::endl;

	DataVector alpha2(numGridPoints);
	alpha2.setAll(0.0);

	DataVector result2(numGridPoints);
	result2.setAll(0.0);

	DataMatrix UpDownMatrix2(numGridPoints, numGridPoints);
	UpDownMatrix2.setAll(0.0);

	myUpDown2 = sg::op_factory::createOperationUpDownTest(*myGrid2);

	std::cout << "start constructing the operator's matrix" << std::endl;
	for (size_t i = 0; i < numGridPoints; i++)
	{
		// init alpha
		alpha2.setAll(0.0);
		alpha2.set(i, 1.0);
		// init result
		result2.setAll(0.0);

		myUpDown2->mult(alpha2, result2);

		// copy data to opartor's matrix
		UpDownMatrix2.setColumn(i, result2);
	}
	std::cout << "finished constructing the operator's matrix" << std::endl;

#ifdef PRINTMATRIX
	std::cout << "The operator's matrix is:" << std::endl << std::endl;
	std::ofstream outfile2;

	outfile2.open("updown.test");

	for (size_t i = 0; i < numGridPoints; i++)
	{
		for (size_t j = 0; j < numGridPoints; j++)
		{
			std::cout << std::scientific << UpDownMatrix2.get(i, j) << " ";
			outfile2 << std::scientific << UpDownMatrix2.get(i, j) << " ";
		}
		std::cout << std::endl;
		outfile2 << std::endl;
	}

	outfile2.close();
#endif

#ifdef PRINTDIFF
	std::string file2 = "XdPhiPhi1D.txt";
	std::ifstream infile2;

	infile2.open(file2.c_str());

	double filedata2 = 0.0;

	for (size_t i = 0; i < numGridPoints; i++)
	{
		for (size_t j = 0; j < numGridPoints; j++)
		{
			infile2 >> filedata2;
			std::cout << (UpDownMatrix2.get(i, j) - filedata2) << " ";
		}
		std::cout << std::endl;
	}

	infile2.close();
#endif

	std::cout << std::endl;

	std::cout << "Test symmetry:" << std::endl;
	double symTest2 = 0.0;
	double tempSym2 = 0.0;
	for (size_t i = 0; i < numGridPoints; i++)
	{
		for (size_t j = 0; j < numGridPoints; j++)
		{
			tempSym2 = fabs(UpDownMatrix2.get(i,j)-UpDownMatrix2.get(j,i));
			if (tempSym2 > symTest2)
			{
				symTest2 = tempSym2;
			}
		}
	}
	std::cout << "Maximum symmetry error: " << symTest2 << std::endl << std::endl;

	std::cout << "Do a Test multiplication:" << std::endl << std::endl;

	alpha2.setAll(0.0);
	alpha2.set(1, 35.0);
	alpha2.set(2, -17.5);
	alpha2.set(4, -7.5);
	alpha2.set(7, -5);
	alpha2.set(14, -1.25);

	result2.setAll(0.0);

	myUpDown2->mult(alpha2, result2);

	std::cout << result2.toString() << std::endl;

	std::cout << std::endl;

	/*sg::IOToolBonnSG* myImporter = new sg::IOToolBonnSG();

		DataVector serAlpha(0);
		std::string serGrid;
		bool hier;

		myImporter->readFile("sparse_bonn.grid", serGrid, serAlpha, hier);

		myImporter->writeFile("sgpp_bonn.grid", *myGrid, alpha, true);

		delete myImporter;*/
	delete myUpDown2;
	delete myGrid2;
#endif
#ifdef WINDOWS
	system("pause");
#endif
	return 0;
}
