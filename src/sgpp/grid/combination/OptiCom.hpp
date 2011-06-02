/*
 * OptiCom.hpp
 *
 *  Created on: Aug 30, 2010
 *      Author: aliz
 */

#ifndef SG_OPTICOM_HPP_
#define SG_OPTICOM_HPP_

#include <iostream>
#include <math.h>
#include <time.h>

#include "sgpp.hpp"

namespace sg {
namespace combigrid {

/** class to contain the methods to calculate the opticom coefficients */
class OptiCom {
	typedef sg::base::GridStorage::index_type::level_type level_t;

public:

  /** the static method to calculate the new coefficients for the combination technique (opticom) */
  static void calc_coefs(sg::base::DataVector &c, FullGridSet* fgs)
  {
	  double A2[4][4]={{2.0,1.0,-2.0,-1.0},{1.0,2.0,-1.0,-2.0},{-2.0,-1.0,2.0,1.0},{-1.0,-2.0,1.0,2.0}};
	  double A1[4][4]={{2.0,-2.0,1.0,-1.0},{-2.0,2.0,-1.0,1.0},{1.0,-1.0,2.0,-2.0},{-1.0,1.0,-2.0,2.0}};
	  int size=fgs->getSize();
	  int dim=2;
	  double* U=new double[size*size];
	  double* E=new double[size];
	  int i,j,k,nr,ii,k0,k1;
	  double sum;
	  sg::base::DataVector coords(dim);
	  vector<level_t>lev(dim);
	  double c1,c2;
	  for (i=0;i<size;i++){
		  for (j=0;j<size;j++)
		  {
			  FullGrid &u=*(fgs->at(i));
			  FullGrid &v=*(fgs->at(j));
			  for (k=0;k<dim;k++)
				  if (u.getLevel()[k]<v.getLevel()[k]) lev[k]=v.getLevel()[k];
				  else lev[k]=u.getLevel()[k];
			  FullGrid* w=FullGrid::createLinearBoundaryFullGrid(dim,&lev);
			//  w->setBoundingBox(u.getBoundingBox());
			  int m=w->getSize();
			  sg::base::DataVector a(m);
			  sg::base::DataVector b(m);
			  for (k=0;k<m;k++)
			  {
				  w->getCoords(k,coords);
				  a[k]=u.eval(coords);
				  b[k]=v.eval(coords);
				  //	cout<<a[k]<<"~="<<b[k]<<",";
			  }
			  int l0=w->length(0);
			  int l1=w->length(1);
			  //cout<<endl;
// 			  if (w->getBoundingBox()!=0){
// 				  c1=1.0/w->getBoundingBox()->getIntervalWidth(0)*(double)(w->length(0)-1)/(w->length(1)-1);
// 				  c2=1.0/w->getBoundingBox()->getIntervalWidth(1)*(double)(w->length(1)-1)/(w->length(0)-1);
// 			  }
// 			  else
// 			  {
				c1=(double)(l0-1)/(l1-1);
				c2=(double)(l1-1)/(l0-1);
// 			  }
			  sum=0;
			  nr=1<<dim;
			  // cout<<l0<<"::"<<l1<<endl;

			  for (k1=0;k1<l1-1;k1++)
			  {
				  for (k0=0;k0<l0-1;k0++){
					  k=k0+k1*l0;
					  for (ii=0; ii < nr; ii++){
						  sum=sum+a[k+ii%2+(ii/2)*l0]*(c1*(A1[ii][0]*b[k]+A1[ii][1]*b[k+1]+A1[ii][2]*b[k+l0]+A1[ii][3]*b[k+l0+1])+c2*(A2[ii][0]*b[k]+A2[ii][1]*b[k+1]+A2[ii][2]*b[k+l0]+A2[ii][3]*b[k+l0+1]));
					  }
					//	cout<<sum<<",";
				  }

			  }
			  U[i*size+j]=sum;//((l0-1)*(l1-1));
			  cout<<U[i*size+j]<<"	";
			  if (i==j) E[i]=sum;
			  delete w;
		  }
		  cout<<"\n";
	  }
	  //    for (int i=0;i<size*size;i++)
	  //    	U[i]=5;

	  //Gaussian elimination
	  for (i=0;i<size-1;i++)
	  {
		  double a=U[i*size+i];
		  for (j=i+1;j<size;j++){
			  double b;
			  if (a!=0)
				  b=U[j*size+i]/a;
			  else
				  b=0.0;
			  U[j*size+i]=0;
			  for (k=i+1;k<size;k++)
				  U[j*size+k]-=b*U[i*size+k];
			  E[j]-=b*E[i];
		  }
	  }

	  for (i=size-1;i>=0;i--)
	  {
		  c[i]=E[i];
		  for (j=size-1;j>i;j--)
			  c[i]-=c[j]*U[i*size+j];
		  if (U[i*size+i]!=0)
			  c[i]/=U[i*size+i];
		  else c[i]=1.0;
		  cout<<"c["<<i<<"]="<<c[i]<<",";
	  }

	  cout<<endl;
	  delete[] U;
	  delete[] E;
	  }
  }; // end of class
} // end namaspace
}

#endif /* SG_OPTICOM_HPP_ */
