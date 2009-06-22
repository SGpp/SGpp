/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#include "tools/classification/ARFFTools.hpp"
#include "exception/file_exception.hpp"
#include <fstream>
#include <stdlib.h>
#include <iostream>

namespace sg
{

ARFFTools::ARFFTools()
{
}

ARFFTools::~ARFFTools()
{
}

size_t ARFFTools::getDimension(std::string tfilename)
{
	std::string line;
	std::ifstream myfile (tfilename.c_str());
	size_t numFound = 0;

	if (myfile.is_open())
	{
		while (! myfile.eof() )
		{
			getline (myfile,line);
			if (line.find("@ATTRIBUTE", 0) != line.npos)
			{
				numFound++;
			}
		}
		myfile.close();
	}
	else
	{
		std::string msg = "Unable to open file: " + tfilename;
		throw new file_exception(msg.c_str());
	}

	// the class is not regarded when getting the dimension
	numFound--;

	return numFound;
}

size_t ARFFTools::getNumberInstances(std::string tfilename)
{
	std::string line;
	std::ifstream myfile (tfilename.c_str());
	size_t numInst = 0;

	if (myfile.is_open())
	{
		while (! myfile.eof() )
		{
			getline (myfile,line);
			if (line.find("@DATA", 0) != line.npos)
			{
				numInst = 0;
			}
			else
			{
				numInst++;
			}
		}
		myfile.close();
	}
	else
	{
		std::string msg = "Unable to open file: " + tfilename;
		throw new file_exception(msg.c_str());
	}

	return numInst;
}

void ARFFTools::readTrainingData(std::string tfilename, DataVector& destination)
{
	std::string line;
	std::ifstream myfile (tfilename.c_str());
	bool data = false;
	size_t instanceNo = 0;

	if (myfile.is_open())
	{
		while (! myfile.eof() )
		{
			getline (myfile,line);
			if (data == true)
			{
				writeNewElement(line, destination, instanceNo);
				instanceNo++;
			}

			if (line.find("@DATA", 0) != line.npos)
			{
				data = true;
			}
		}
		myfile.close();
	}
	else
	{
		std::string msg = "Unable to open file: " + tfilename;
		throw new file_exception(msg.c_str());
	}
}

void ARFFTools::readClasses(std::string tfilename, DataVector& destination)
{
	std::string line;
	std::ifstream myfile (tfilename.c_str());
	bool data = false;
	size_t instanceNo = 0;

	if (myfile.is_open())
	{
		while (! myfile.eof() )
		{
			getline (myfile,line);
			if (data == true)
			{
				writeNewClass(line, destination, instanceNo);
				instanceNo++;
			}

			if (line.find("@DATA", 0) != line.npos)
			{
				data = true;
			}
		}
		myfile.close();
	}
	else
	{
		std::string msg = "Unable to open file: " + tfilename;
		throw new file_exception(msg.c_str());
	}
}

void ARFFTools::writeNewElement(std::string& instance, DataVector& destination, size_t instanceNo)
{
	size_t cur_pos = 0;
	size_t cur_find = 0;
	size_t dim = destination.getDim();
	std::string cur_value;
	double dbl_cur_value;

	for (size_t i = 0; i < dim; i++)
	{
		cur_find = instance.find(",", cur_pos);
		cur_value = instance.substr(cur_pos, cur_find-cur_pos);
		dbl_cur_value = atof(cur_value.c_str());
		destination.set((instanceNo*dim) + i, dbl_cur_value);
		cur_pos = cur_find + 1;
	}
}

void ARFFTools::writeNewClass(std::string& instance, DataVector& destination, size_t instanceNo)
{
	size_t cur_pos = instance.find_last_of(",");
	std::string cur_value = instance.substr(cur_pos+1);
	double dbl_cur_value = atof(cur_value.c_str());
	destination.set(instanceNo, dbl_cur_value);
}

void ARFFTools::writeAlpha(std::string tfilename, DataVector& source)
{

}

void ARFFTools::readAlpha(std::string tfilename, DataVector& destination)
{

}

}
