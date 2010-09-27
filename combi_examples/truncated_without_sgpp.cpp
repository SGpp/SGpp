/*
 * clasical_combi.cpp
 *
 *  Created on: Aug 20, 2010
 *      Author: Aliz Nagy
 */

#include <iostream>
#include <math.h>
#include "data/DataVector.hpp"
#include "grid/combination/FullGridSet.hpp"
#include "grid/combination/FullGridSet.cpp"
using namespace sg;
const double PI = 3.14159265;
/**
 * The function we will try to interpolate
 * */
double f(double x0, double x1, double x2)
{
	return 1.0+(0.25*(x0-0.7)*(x0-0.7)+2.0)+(0.25*(x1-0.7)*(x1-0.7)+2.0)+(0.25*(x2-0.7)*(x2-0.7)+2.0);
}
int main()
{
	int dim=3;
	int level=4;
    int l_user=2;
    //Create the set of fullgrids into which the grid decomposes
    FullGridSet fgs(dim,level,l_user);
    cout<<"Number of full grids:"<<fgs.getSize()<<endl;
    cout<<"The fullgrids:"<<endl;
    /**
     * This prints the levels of all fullgrids
     */
    fgs.printGridSet();
    //Create a new datavector which contains the coordinates of the point we want to interpolate
    DataVector p(dim);
    p[0] = 0.91;
    p[1] = 0.23;
    p[2] = 0.76;
    /**
     * Fill the full grids with the function values
     * */
    for (int i=0;i<fgs.getSize();i++)
           {
			size_t m=fgs[i].getSize();
           	FullGrid *fg=fgs.at(i);//works also with &(fgs[i])
           	for (int j=0;j<m;j++)
           	       	{
					 (*fg)[j]=f(fg->getCoord(0,j),fg->getCoord(1,j),fg->getCoord(2,j));
           	       	}
           	/**
           	 * *Evaluates the fullgrid in an arbitrary point, and assigns the resulting value to the field variable value of the fullgrid,
           	 *  The same value is returned by the function and can be accesed later through the function call fullgrid.val()
             */
			fg->eval(p);
           }
    /**
     * Combines the interpolated results on the fullgrids into one value, which in case of function interpolation equals the value on the sparse grid itself
     * We print the interpolation value as Uc(p)
     */
    cout<<"Uc("<<p[0]<<","<<p[1]<<","<<p[2]<<")="<<fgs.combinedResult()<<endl;
    /**
     * We print the real value of the function in the given point
     * */
    cout<<"f("<<p[0]<<","<<p[1]<<","<<p[2]<<")="<<f(p[0],p[1],p[2])<<endl;
    return 0;
}
