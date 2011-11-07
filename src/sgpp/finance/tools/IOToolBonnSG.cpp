/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "finance/tools/IOToolBonnSG.hpp"
#include "base/exception/tool_exception.hpp"

#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <list>

namespace sg
{
namespace pde
{

IOToolBonnSG::IOToolBonnSG()
{
}

IOToolBonnSG::~IOToolBonnSG()
{
}

void IOToolBonnSG::readFile(std::string tFilename, std::string& sgppSerialization, sg::base::DataVector& alpha, bool& ishierarchized)
{
	size_t dim;
	size_t gridpoints;
	std::string curLine;
	std::string GridPointSerialization;
	char* curLine_cstr;
	char* token;

	std::ifstream file;
	char line[8192];

	file.open(tFilename.c_str());
	if(!file.is_open())
	{
		throw new sg::base::tool_exception("IOToolBonnSG::readFile : The specified grid file doesn't exists!");
	}

	// parse the global grid information
	file.getline(line,4096);
	curLine.assign(line);

	curLine_cstr = (char *)curLine.c_str();
	token = strtok(curLine_cstr, ",");

	// read dimension
	char* strDim_cstr;
	strDim_cstr = strtok(0, ",");

	// read number grid points
	char* strSize_cstr;
	strSize_cstr = strtok(0, ",");

	// read coefficient type
	char* strType_cstr;
	strType_cstr = strtok(0, ",");

	token = strtok(strDim_cstr, "=");
	token = strtok(0, "=");
	dim = atoi(token);

	token = strtok(strSize_cstr, "=");
	token = strtok(0, "=");
	gridpoints = atoi(token);

	token = strtok(strType_cstr, "=");
	token = strtok(0, "=");

	std::string strType;
	strType.assign(token);

	if (strType == "BASIS_NODAL")
	{
		ishierarchized = false;
	}
	else if (strType == "BASIS_HAT_HIER_LINEAR")
	{
		ishierarchized = true;
	}
	else
	{
		throw new sg::base::tool_exception("IOToolBonnSG::readFile : The specified grid file contains an unknown basis!");
	}

	sgppSerialization = "";
	GridPointSerialization = "";
	std::stringstream gridinfo;
	gridinfo << "linearTrapezoidBoundary" << std::endl;
	gridinfo << "1" << std::endl;
	gridinfo << "4 " << dim << " " << gridpoints;
	sgppSerialization.append(gridinfo.str());

	alpha.resize(gridpoints);

	std::list<double>* coords = new std::list<double>[dim];		// Array of Lists that store all coordiantes
	std::string curCoords;				// contains coordinates of current grid point

	// parse the grid points and dertermine bounding box
	for (size_t i = 0; i < gridpoints; i++)
	{
		// read next line from file
		file.getline(line,8192);

		curLine.assign(line);

		std::string SGPoint;
		std::string SGLevel;
		std::string SGIndex;
		std::string curAlpha;
		std::stringstream streamSGPoint;
		std::stringstream streamLevel;
		std::stringstream streamIndex;

		char* curLine_cstr = (char *)curLine.c_str();
		token = strtok(curLine_cstr, "=");

		// contains (l,l,...,l|i,i,...,i)
		SGPoint.assign(token);

		token = strtok (0, "=");

		// read cooradinates
		curCoords.assign(token);

		token = strtok (0, "=");
		token = strtok (0, "=");

		// contains number
		curAlpha.assign(token);

		// remave brackets
		SGPoint = SGPoint.substr(SGPoint.find("(")+1, (SGPoint.find(")") - SGPoint.find("(")) - 1);

		// Split level and index
		char* SGPoint_cstr;
		char* SGLevel_cstr;
		char* SGIndex_cstr;
		SGPoint_cstr = (char *)SGPoint.c_str();

		// contains all levels of sg::base::Grid Point
		SGLevel_cstr = strtok(SGPoint_cstr, "|");
		// conatins all indeced of sg::base::Grid Point
		SGIndex_cstr = strtok (0, "|");

		size_t* levels = new size_t[dim];
		size_t* indeces = new size_t[dim];

		token = strtok(SGLevel_cstr, ",");
		levels[0] = atoi(token);
		for (size_t j = 1; j < dim; j++)
		{
			token = strtok(0, ",");
			levels[j] = atoi(token);
		}
		token = strtok(SGIndex_cstr, ",");
		indeces[0] = atoi(token);
		for (size_t j = 1; j < dim; j++)
		{
			token = strtok(0, ",");
			indeces[j] = atoi(token);
		}

		// reorder level and index and build SGpp's index definition
		streamSGPoint << std::endl << dim << std::endl;
		for (size_t j = 0; j < dim; j++)
		{
			streamSGPoint << levels[j] << " " << indeces[j] << " ";
		}

		delete[] levels;
		delete[] indeces;

		// append definition to string, write alpha entry
		GridPointSerialization.append(streamSGPoint.str());
		alpha.set(i, atof(curAlpha.c_str()));

		// Determine Coordinates
		// remove brackets
		curCoords = curCoords.substr(curCoords.find("(")+1 ,  (curCoords.find(")") - curCoords.find("(")) - 1);
		// Split Coordinates
		char* curCoords_cstr = (char *)curCoords.c_str();
		token = strtok(curCoords_cstr, ",");
		coords[0].push_back(atof(token));
		for (size_t j = 1; j < dim; j++)
		{
			token = strtok(0, ",");
			coords[j].push_back(atof(token));
		}
	}

	// Determin Bounding Box
	std::stringstream streamBoudingBox;
	typedef std::list<double>::iterator iter;

	streamBoudingBox << std::endl;

	for (size_t j = 0; j < dim; j++)
	{
		coords[j].sort();

		size_t counter = 0;

		for (iter i = coords[j].begin(); i != coords[j].end(); i++)
		{
			if (counter == 0 || counter == (gridpoints - 1))
			{
				streamBoudingBox << std::scientific << *i << " ";
			}
			counter++;
		}
		streamBoudingBox << "0 0 ";
	}

	sgppSerialization.append(streamBoudingBox.str());
	sgppSerialization.append(GridPointSerialization);

	//std::cout << sgppSerialization << std::endl;
	//std::cout << alpha.toString() << std::endl;
}

void IOToolBonnSG::writeFile(std::string tFilename, sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, bool ishierarchized)
{
	std::ofstream fout;

	fout.open(tFilename.c_str());
	if(!fout.is_open())
	{
		throw new sg::base::tool_exception("IOToolBonnSG::writeFile : The specified grid file could not be open for writing!");
	}

	// Write global gird information
	fout << "datastructure=hash::adp,dim=" << SparseGrid.getStorage()->dim() <<  ",basissize=" << SparseGrid.getStorage()->size() << ",basis=";

	if (ishierarchized == true)
	{
		fout << "BASIS_HAT_HIER_LINEAR";
	}
	else
	{
		fout << "BASIS_NODAL";
	}

	// serialize the gridpoints
	for (size_t i = 0; i < SparseGrid.getStorage()->size(); i++)
	{
		std::stringstream strlevel;
		std::stringstream strindex;
		sg::base::GridStorage::index_type::level_type level;
		sg::base::GridStorage::index_type::index_type index;
		bool isBoundary = false;

		// print the sparse grid index
		for (size_t j = 0; j <  SparseGrid.getStorage()->dim(); j++)
		{
			SparseGrid.getStorage()->get(i)->get(j, level, index);
			if (level == 0)
			{
				isBoundary = true;
			}
			strlevel << level;
			strindex << index;

			if (j < (SparseGrid.getStorage()->dim()-1))
			{
				strlevel << ",";
				strindex << ",";
			}
		}
		fout << std::endl << "(" << strlevel.str() << "|" << strindex.str() << ")	= ";

		// print the 'real' coordinates
		fout << "(";
		for (size_t j = 0; j <  SparseGrid.getStorage()->dim(); j++)
		{
			fout << std::scientific << SparseGrid.getStorage()->get(i)->getCoordBB(j, SparseGrid.getBoundingBox()->getIntervalWidth(j), SparseGrid.getBoundingBox()->getIntervalOffset(j));

			if (j < (SparseGrid.getStorage()->dim()-1))
			{
				fout << ",";
			}
		}
		fout << ")	= ";

		// print point description
		std::stringstream hash;  // use this to avoid scientific format
		hash << i;
		if (isBoundary == true)
		{
			fout << "Boundary:";
		}
		else
		{
			fout << "Inner:";
		}
		fout << hash.str() << "	=";

		// print alpha
		fout << std::scientific << alpha.get(i);
	}

	fout.close();
}

}
}
