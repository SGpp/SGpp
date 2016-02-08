///* ****************************************************************************
//* Copyright (C) 2011 Technische Universitaet Muenchen                         *
//* This file is part of the SG++ project. For conditions of distribution and   *
//* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
//**************************************************************************** */
//// @author Janos Benk (benk@in.tum.de)
//// modified by: PeTz
//
//#include "combigrid/plotter/GridPlotter.hpp"
//
//using namespace combigrid;
//
//template
//
//static void GridPlotter<_FGTp,_CFTp>::plotFullGrid(const std::string& filePath, const FullGrid<_FGTp>* fg,
//		std::vector<double>& globalCoord_in, int resolution = 0) {
//
//	combigrid::Evaluable<_FGTp, _CFTp> obj(fg);
//	int dim = fg->getDimension();
//	plotObject(dim, filePath, &obj, fg->getDomain(), globalCoord_in,
//			resolution);
//}
//
///** plot one combination grid */
//template<typename _FGTp,template _CFTp>
//static void GridPlotter<_FGTp,_CFTp>::plotCombiGrid(const std::string& filePath,
//		const CombiGrid_demo<_FGTp, _CFTp>* cg,
//		std::vector<double>& globalCoord_in, int resolution = 0) {
//
//	combigrid::Evaluable<_FGTp, _CFTp> obj(cg);
//	int dim = cg->getFullGrid(0)->getDimension();
//	plotObject(dim, filePath, &obj, cg->getDomain(), globalCoord_in,
//			resolution);
//
//}
//
//template<typename _FGTp,template _CFTp>
//static void GridPlotter<_FGTp,_CFTp>::plotObject(int dim, const std::string& filePath,
//		const combigrid::Evaluable<_FGTp, _CFTp>* obj, const GridDomain* domain,
//		std::vector<double>& globalCoord_in, int resolution) {
//	std::vector<_CFTp> result(0);
//	std::vector<double> globalCoord = globalCoord_in;
//	bool isScaled = false;
//
//	if (domain != 0) {
//		isScaled = domain->get1DDomain(0).isAxisScaled();
//	}
//
//	if (((domain == 0) || (resolution > 0)) || (!isScaled)) {
//
//		if (dim == 1) {
//			// if the domain is not existent then just use a regular resolution
//			if (resolution <= 0)
//			resolution = 50;
//
//			double minX = 0.0;
//			double maxX = 1.0;
//
//			if (domain != 0) {
//				minX = domain->get1DDomain(0).getMinDomain();
//				maxX = domain->get1DDomain(0).getMaxDomain();
//			}
//
//			result.resize(resolution);
//
//			// loop and evaluate points
//			for (int ii = 0; ii < resolution; ii++) {
//				globalCoord[0] = minX
//				+ (maxX - minX) * (double(ii) / double(resolution - 1));
//				result[ii] = obj->eval(globalCoord);
//			}
//
//			// writing file
//			ofstream myfile;
//			myfile.open(filePath.c_str());
//			myfile << "X = [ " << minX;
//
//			for (int ii = 1; ii < resolution; ii++) {
//				myfile << " , "
//				<< (minX
//						+ (maxX - minX)
//						* (double(ii) / double(resolution - 1)));
//			}
//
//			myfile << "]; \n ";
//			myfile << "res = [ " << result[0];
//
//			for (int ii = 1; ii < resolution; ii++) {
//				myfile << " , " << result[ii];
//			}
//
//			myfile << "]; \n ";
//			myfile << " plot(X,res); \n " << result[0];
//			myfile.close();
//		} else {
//			double minX = 0.0;
//			double maxX = 1.0;
//			double minY = 0.0;
//			double maxY = 1.0;
//
//			if (domain != 0) {
//				minX = domain->get1DDomain(0).getMinDomain();
//				maxX = domain->get1DDomain(0).getMaxDomain();
//				minY = domain->get1DDomain(1).getMinDomain();
//				maxY = domain->get1DDomain(1).getMaxDomain();
//			}
//
//			result.resize(resolution * resolution);
//
//			// loop and evaluate points
//			for (int ii = 0; ii < resolution; ii++) {
//				for (int jj = 0; jj < resolution; jj++) {
//					globalCoord[0] = minX
//					+ (maxX - minX)
//					* (double(ii) / double(resolution - 1));
//					globalCoord[1] = minY
//					+ (maxY - minY)
//					* (double(jj) / double(resolution - 1));
//					result[ii * resolution + jj] = obj->eval(globalCoord);
//				}
//			}
//
//			ofstream myfile;
//			myfile.open(filePath.c_str());
//			myfile << "X = [ " << minX;
//
//			for (int ii = 1; ii < resolution; ii++) {
//				myfile << " , "
//				<< (minX
//						+ (maxX - minX)
//						* (double(ii) / double(resolution - 1)));
//			}
//
//			myfile << "]; \n ";
//			myfile << "Y = [ " << minY;
//
//			for (int ii = 1; ii < resolution; ii++) {
//				myfile << " , "
//				<< (minY
//						+ (maxY - minY)
//						* (double(ii) / double(resolution - 1)));
//			}
//
//			myfile << "]; \n ";
//			myfile << "res = [ " << result[0];
//
//			for (int ii = 1; ii < resolution * resolution; ii++) {
//				if ((ii % resolution) == 0) {
//					myfile << " ; " << result[ii];
//				} else {
//					myfile << " , " << result[ii];
//				}
//			}
//
//			myfile << "]; \n ";
//			myfile << "[x,y]=meshgrid(Y,X);\n";
//			myfile << " surf(x,y,res); \n ";
//			myfile.close();
//		}
//	} else {
//		if (dim == 1) {
//			result.resize(domain->get1DDomain(0).axisScaling().size());
//
//			// loop and evaluate points
//			for (unsigned int ii = 0;
//					ii < domain->get1DDomain(0).axisScaling().size(); ii++) {
//				globalCoord[0] = domain->get1DDomain(0).axisScaling()[ii];
//				result[ii] = obj->eval(globalCoord);
//			}
//
//			// writing file
//			ofstream myfile;
//			myfile.open(filePath.c_str());
//			myfile << "X = [ " << domain->get1DDomain(0).axisScaling()[0];
//
//			for (unsigned int ii = 1;
//					ii < domain->get1DDomain(0).axisScaling().size(); ii++) {
//				myfile << " , " << domain->get1DDomain(0).axisScaling()[ii];
//			}
//
//			myfile << "]; \n ";
//			myfile << "res = [ " << result[0];
//
//			for (unsigned int ii = 1;
//					ii < domain->get1DDomain(0).axisScaling().size(); ii++) {
//				myfile << " , " << result[ii];
//			}
//
//			myfile << "]; \n ";
//			myfile << " plot(X,res); \n ";
//			myfile.close();
//		} else {
//			result.resize(
//					domain->get1DDomain(0).axisScaling().size()
//					* domain->get1DDomain(1).axisScaling().size());
//			// loop and evaluate points
//			double res = 0.0;
//
//			for (unsigned int ii = 0;
//					ii < domain->get1DDomain(0).axisScaling().size(); ii++) {
//				for (unsigned int jj = 0;
//						jj < domain->get1DDomain(1).axisScaling().size();
//						jj++) {
//					globalCoord[0] = domain->get1DDomain(0).axisScaling()[ii];
//					globalCoord[1] = domain->get1DDomain(1).axisScaling()[jj];
//					res = obj->eval(globalCoord);
//					result[ii * domain->get1DDomain(1).axisScaling().size() + jj] =
//					res;
//				}
//			}
//
//			ofstream myfile;
//			myfile.open(filePath.c_str());
//			myfile << "X = [ " << domain->get1DDomain(0).axisScaling()[0];
//
//			for (unsigned int ii = 1;
//					ii < domain->get1DDomain(0).axisScaling().size(); ii++) {
//				myfile << " , " << domain->get1DDomain(0).axisScaling()[ii];
//			}
//
//			myfile << "]; \n ";
//			myfile << "Y = [ " << domain->get1DDomain(1).axisScaling()[0];
//
//			for (unsigned int ii = 1;
//					ii < domain->get1DDomain(1).axisScaling().size(); ii++) {
//				myfile << " , " << domain->get1DDomain(1).axisScaling()[ii];
//			}
//
//			myfile << "]; \n ";
//			myfile << "res = [ " << result[0];
//
//			for (unsigned int ii = 1;
//					ii
//					< domain->get1DDomain(0).axisScaling().size()
//					* domain->get1DDomain(1).axisScaling().size();
//					ii++) {
//				if ((ii % domain->get1DDomain(1).axisScaling().size()) == 0) {
//					myfile << " ; " << result[ii];
//				} else {
//					myfile << " , " << result[ii];
//				}
//			}
//
//			myfile << "]; \n ";
//			myfile << "[x,y]=meshgrid(Y,X);\n";
//			myfile << " surf(x,y,res); \n ";
//			myfile.close();
//		}
//	}
//}
//
///**
//* In order to include different implementations of the Evaluable template class
//* into the compiled combigrid static/shared libraries, it is necessary to append the following template
//* class declarations in this  cpp file:
//*
//* 	"template class combigrid::Evaluable<type1,type2>;"
//*
//* where type1 specifies the template type for the fullgrid elements and type2 is the template type
//* for the combigrid coefficients!
//*
//* IMPORTANT: For logical consistency and correct program execution, type1 and type2,
//* i.e. the fullgrid elements type, and the coefficients type should satisfy an imposed by us
//* "type precedence relation". That is, the fullgrid(type1) template type should have LOWER or EQUAL precedence
//* to the coefficients template type(type2)! This is necessary for correct function evaluation on the fullgrids.
//*
//* In somewhat more formal notation we require that:
//* 	precedence(type2) >= precedence(type1)
//*
//* 	If this condition is violated, during class initialization a runtime exception is thrown and the
//* 	execution of the program is aborted! (see combigrid_utils.hpp -> check_type_precedence() function!
//*
//* In this context the following class declarations should be valid:
//*  combigrid::Evaluable<float,float>
//*  combigrid::Evaluable<float,double>;
//* end the following should not be valid
//*	combigrid::Evaluable<double,float>;
//*	combigrid::Evaluable<float,int>;
//*
//*
//* */
//
//// FULLGRID type is float, coeffs type is also float
//template class combigrid::Evaluable<float, float>;
//// FULLGRID type is float coeffs type is double
//template class combigrid::Evaluable<float, double>;
//// FULLGRID type is double, coefs type is also double!
//template class combigrid::Evaluable<double, double>;
//
//// add more declarations at your hearth's will !!!!
//
//template class combigrid::GridPlotter<float, float>;
//// FULLGRID type is float coeffs type is double
//template class combigrid::GridPlotter<float, double>;
//// FULLGRID type is double, coefs type is also double!
//template class combigrid::GridPlotter<double, double>;
//// add more declarations at your hearth's will !!!!
