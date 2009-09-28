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

#include "application/common/ScreenOutput.hpp"

namespace sg
{

#ifdef WINDOWS
ScreenOutput::ScreenOutput()
{
	hCon = GetStdHandle(STD_OUTPUT_HANDLE);
	first_run = true;
}

void ScreenOutput::update(int progress, std::string status)
{
	int i;
	GetConsoleScreenBufferInfo(hCon,&info);

	if(!first_run)
	{
		pos.Y = info.dwCursorPosition.Y - 3;
		pos.X = 0;
		SetConsoleCursorPosition(hCon,pos);
	}
	else
	{
		first_run = false;
	}

	std::cout << "[";

	for(i=0; i<(progress/2); i++)
	{
		std::cout << '\xdb';
	}
	for(; i<50; i++)
	{
		std::cout << " ";
	}

	std::cout << "]  " << progress << "%" << std::endl << std::endl;
	std::cout << status << "       " << std::endl;
}
#endif

#ifndef WINDOWS
ScreenOutput::ScreenOutput()
{
	first_run = true;
}

void ScreenOutput::update(size_t progress, std::string status)
{
	size_t i;

	if(!first_run)
	{
		std::cout << "\033[4A" << std::endl;
	}
	else
	{
		first_run = false;
	}

	std::cout << "[";

	for(i=0; i<(progress/2); i++)
	{
		std::cout << "=";
	}
	for(; i<50; i++)
	{
		std::cout << " ";
	}

	std::cout << "]  " << progress << "%" << std::endl << std::endl;
	std::cout << status << "          " << std::endl;
}
#endif

ScreenOutput::~ScreenOutput()
{
}

void ScreenOutput::writeTitle(std::string appTitle, std::string appAuthor)
{
	std::cout << std::endl << std::endl;

	std::string empty = "                                                      ";
	appTitle.append(empty);
	appAuthor.append(empty);

	appTitle = appTitle.substr(0, empty.length()-2);
	appAuthor = appAuthor.substr(0, empty.length());

	std::cout << appTitle << "########  ##########" << std::endl;
	std::cout << "=========================================           " << "  ##  ##  ##  ##  ##" << std::endl;
	std::cout << "                                                      ##  ##  ##  ##  ##" << std::endl;
	std::cout << appAuthor << "##  ##  ##  ##  ##" << std::endl;
	std::cout << "                                                      ##  ######  ##  ##" << std::endl << std::endl;
}

void ScreenOutput::writeHelp(std::string helpText)
{
	std::cout << helpText << std::endl << std::endl;
}

void ScreenOutput::writeStartSolve(std::string text)
{
	std::string line;
	size_t line_length = (size_t)text.length();
	while(line_length > 0)
	{
		line.append("-");
		line_length--;
	}

	std::cout << std::endl;
	std::cout << text << std::endl;
	std::cout << line << std::endl;
	std::cout << std::endl;
}

void ScreenOutput::writeEmptyLines(size_t numLines)
{
	for (size_t i = 0; i < numLines; i++)
	{
		std::cout << std::endl;
	}
}

}
