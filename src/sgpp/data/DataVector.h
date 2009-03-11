/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2007 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/


#ifndef DATAVECTOR_H_
#define DATAVECTOR_H_

#include <string>
#include <vector>

/**
 * a data holding class.
 */
class DataVector
{
public:
	DataVector(int size);
	DataVector(int size, int dim);
	DataVector(DataVector& vec);
	DataVector(double* input, int size, int dim);

	void resize(int size);
	void addSize(int add);
	int addValue();

	void setAll(double value);

	void copyFrom(const DataVector& vec);
	void copySmall(const DataVector& vec);
	DataVector& operator=(const DataVector& vec);
	inline double& operator[](int i)
	{
		return data[i];
	};

	double get(int i) const;
	void set(int i, double value);

 	void getRow(int row, DataVector& vec);
 	void setRow(int row, DataVector& vec);
 	void getColumn(int col, DataVector& vec);
 	void setColumn(int col, DataVector& vec);

	void add(DataVector& vec);
	void sub(DataVector& vec);
	void mult(double scalar);

	void sqr();
	double sum();

	void axpy(double alpha, DataVector& x);

	void getLine(int row, DataVector& vec);
	void getLine(int row, std::vector<double>& vec);

	double dotProduct(DataVector& vec);

	int getSize();
	int getDim();
	int getTotalSize();
	inline int getUnused() { return unused; };

	void partitionClasses(double border);
	void normalizeDimension(int d);
	void normalizeDimension(int d, double border);

	double min(int d);
	double max(int d);
	void minmax(int d, double* min, double* max);

	void toString(std::string& text);

	double* getPointer();

	virtual ~DataVector();

private:
	DataVector();

	/// array to store the data
	double* data;
	/// number of elements in the data vector
	int size;
	/// number of dimensions of one element in this vector
	int dim;
	/// unused slots in this data vector
	int unused;
};

#endif /*DATAVECTOR_H_*/
