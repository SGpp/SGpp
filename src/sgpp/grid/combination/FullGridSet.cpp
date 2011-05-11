#include "FullGridSet.hpp"

#include <string.h>
#include <stdlib.h>
using namespace sg::base;

namespace sg{
namespace combigrid {

void FullGridSet::generate(size_t dim,size_t n)
{
        gridsWithBonudaryPoints_ = true;
		vector<level_t> v(0);
		size=0;
		int combi,oldsize;
        ///generates fullgrids with norm1 up to level to get all grid points
		for (size_t q=0;q<dim;q++)
		{							
	    	combi=combination(dim-1,q);
	    	oldsize=size;
	     	if (q%2!=0) combi=-combi;
			getsums(&v,dim,n-q);
			for (size_t i=oldsize;i<size;i++)
				coefs.push_back(combi);
			
		}          
		gridType=0;
}

void FullGridSet::generate(size_t dim,size_t n,const char *type)
{
		if (strcmp(type,"linearBoundary")==0){
		    gridsWithBonudaryPoints_ = true;
			generate(dim,n);
		}
        else
        {
        	if (strstr(type,"linear")!=0) {
        		vector<level_t> v(0);
        		size=0;
        		int combi,oldsize;
        		if (strstr(type,"Trapezoid")!=0){
            		 gridsWithBonudaryPoints_ = true;
                     n=n+dim-1;///generate fullgrids with norm1 up to level+dim-1 to get all gridpoints	                                
                     for (size_t q=0;q<dim;q++)
                     {
                    	if ((n-q)>=dim){
                    	     combi=combination(dim-1,q);
                    	     oldsize=size;
                    	  	if (q%2!=0) combi=-combi;
                    	   	 getTrapezoidsums(&v,dim,n-q,1);
                        	 for (size_t i=oldsize;i<size;i++)
                        	 		coefs.push_back(combi);
                         }
                     }
                     gridType=2;          
        		}
        		else{
        			 gridsWithBonudaryPoints_ = false;
                     n=n+dim-1;///generate fullgrids with norm1 up to level+dim-1 to get all gridpoints
                     for (size_t q=0;q<dim;q++)
                     {
                            if ((n-q)>=dim){
                            	combi=combination(dim-1,q);
                             	if (q%2!=0) combi=-combi;
                            	oldsize=size;
                            	getInnersums(&v,dim,n-q);
                            	for (size_t i=oldsize;i<size;i++)
                            	      coefs.push_back(combi);
                            }
                     }
                    gridType=1;              
                 }		  		
			}
        	else if (strcmp(type,"squareRoot")==0){
    	         vector<level_t> v(0);
    	         size=0;
        	     getSquare(&v,dim,n);
        	     gridType=3;
        	}
        	else cout<<"No combination technique implemented for this type of grid";
        }
}

void FullGridSet::generate(size_t dim,size_t n,size_t l_user)
{
	  vector<level_t> v(0);
	  size=0;
	  int combi,oldsize;
	  n=n+(dim-1)*l_user;///generate fullgrids with norm1 up to level+dim-1 to get all gridpoints
	  for (size_t q=0;q<dim;q++)
	  {
		 if ((n-q)>=dim*l_user){
			 combi=combination(dim-1,q);
			 if (q%2!=0) combi=-combi;
			 oldsize=size;
			 getTrapezoidsums(&v,dim,n-q,l_user);
			 for (size_t i=oldsize;i<size;i++)
			         coefs.push_back(combi);
		 }
	  }
	  gridType=2;
}

void FullGridSet::generate(size_t dim,size_t n,size_t *l_user)
{
	  vector<level_t> v(0);
	  size=0;
	  size_t max=l_user[0];
	  size_t sum=max;
	  int combi,oldsize;
	  for (size_t i=1;i<dim;i++)
	  {
			  sum+=l_user[i];
			  if (l_user[i]>max) max=l_user[i];
	  }

	  n=n+sum-max;
	  for (size_t q=0;q<dim;q++)
	  {
		 if ((n-q)>=sum){
			 combi=combination(dim-1,q);
			 if (q%2!=0) combi=-combi;
			 oldsize=size;
			 getTrapezoidsums(&v,dim,n-q,l_user);
			 for (size_t i=oldsize;i<size;i++)
			           coefs.push_back(combi);
		 }
	  }
	  gridType=2;
}

void FullGridSet::generate(size_t dim, size_t* levels, const char *type,size_t l_user)
{
	  vector<level_t>* v=new vector<level_t>(0);
	  size=0;
		if (strcmp(type,"squareRoot")==0)
		{
			for (size_t i=0;i<dim;i++)
			{
					  v->clear();
					  for (size_t j=0;j<dim;j++)
					  {
						  if (i!=j) v->push_back(levels[j]/2);
						  else v->push_back(levels[j]);
					  }
					  if (gridsWithBonudaryPoints_ )
						{grids.push_back(FullGrid::createLinearBoundaryFullGrid(v->size(),v));}
					  else
						{grids.push_back(FullGrid::createLinearFullGrid(v->size(),v));}
					  coefs.push_back(1);
			}
			v->clear();
			for (size_t j=0;j<dim;j++)
			{
					  v->push_back(levels[j]/2);
			}
			if (gridsWithBonudaryPoints_ )
				{grids.push_back(FullGrid::createLinearBoundaryFullGrid(v->size(),v));}
			else
				{grids.push_back(FullGrid::createLinearFullGrid(v->size(),v));}
			coefs.push_back(1.0-dim);
			size=dim+1;
		}
		else if (strstr(type,"linear")!=0)
		{
			size_t n=levels[0];
			size_t i,q;
			int combi,oldsize;
			double* ratio=new double[dim];
			for (i=1;i<dim;i++){
				if (levels[i]>n)
					n=levels[i];
			}
			for (i=0;i<dim;i++){
				ratio[i]=(double)n/levels[i];
			}
			if (strcmp(type,"linearBoundary")==0) {
				        gridsWithBonudaryPoints_ = true;
						for (q=0;q<dim;q++)
						{
					    	combi=combination(dim-1,q);
					    	oldsize=size;
					     	if (q%2!=0) combi=-combi;
							getTrapezoidsums(v,dim,n-q,0,ratio);
							for (i=oldsize;i<size;i++)
								coefs.push_back(combi);
						}
						gridType=0;
						removeDuplicates();
			}
			else if (strcmp(type,"linearTrapezoidBoundary")==0){
				         gridsWithBonudaryPoints_ = true;
					     n=n-1;
						 for (q=0;q<dim;q++)
							 {
									combi=combination(dim-1,q);
									oldsize=size;
									if (q%2!=0) combi=-combi;
									 getTrapezoidsums(v,dim,n-q,1,ratio);
									 for (i=oldsize;i<size;i++)
											coefs.push_back(combi);
							  }
						 gridType=2;
						 removeDuplicates();
			}
			else if (strcmp(type,"modlinearTrapezoidBoundary")==0){
				      gridsWithBonudaryPoints_ = true;
					  n=n-l_user;
					  for (q=0;q<dim;q++)
					  {

							 combi=combination(dim-1,q);
							 if (q%2!=0) combi=-combi;
							 oldsize=size;
							 getTrapezoidsums(v,dim,n-q,l_user,ratio);
							 for (i=oldsize;i<size;i++)
							                coefs.push_back(combi);

					  }
					  gridType=2;
					  removeDuplicates();
			}
			else {
				 gridsWithBonudaryPoints_ = false;
				 n=n-1;
				 for (q=0;q<dim;q++)
				 {
						combi=combination(dim-1,q);
						if (q%2!=0) combi=-combi;
						oldsize=size;
						getInnersums(v,dim,n-q,ratio);
						for (i=oldsize;i<size;i++)
									   coefs.push_back(combi);
				 }
				 gridType=1;
			}
		}
		else cout<<"Grid type not supported";
}

void FullGridSet::generate(size_t dim, size_t* levels, size_t* l_user)
{
		size_t n=levels[0];
		size_t i,q;
		size_t max=l_user[0];
		vector<level_t> v(0);

		//std::cout << " FullGridSet::generate boundary:" << gridsWithBonudaryPoints_ << std::endl;

		size=0;
		int combi,oldsize;
		double* ratio=new double[dim];
		for (i=1;i<dim;i++){
			if (levels[i]>n)
				n=levels[i];
			if (l_user[i]>max)
				max=l_user[i];
		}
		for (i=0;i<dim;i++){
			ratio[i]=(double)n/levels[i];
		}
		n=n-max;
		for (q=0;q<dim;q++)
		{
				combi=combination(dim-1,q);
				oldsize=size;
				if (q%2!=0) combi=-combi;
				getTrapezoidsums(&v,dim,n-q,l_user,ratio);
				for (i=oldsize;i<size;i++)
						coefs.push_back(combi);
		 }
		 gridType=2;
		 removeDuplicates();
}

void FullGridSet::removeDuplicates()
{
	size_t i=0;
	size_t j;
	vector<FullGrid*>::iterator gt=grids.begin();
	vector<double>::iterator ct=coefs.begin();
	while(i<size-1){
		j=i+1;
		while(j<size)
		{
			if (grids.at(i)->equals(*(grids.at(j)))){
					grids.erase(gt+j);
					coefs[i]+=coefs[j];
					coefs.erase(ct+j);
					size--;
			}
			else j++;
		}
		if (coefs[i]==0){
			grids.erase(gt+i);
			coefs.erase(ct+i);
			size--;
			if (i==0){
				gt=grids.begin();
				ct=coefs.begin();
			}
		}
		else i++;
	}
}

void FullGridSet::getsums(vector<level_t> *v,size_t dim,size_t sum)
{
/* Takes recursively every possible combination of numbers which add up to sum creating a linear boundary grid for each one */
		if (dim==1)
		{
			    v->push_back(sum);
			    if (gridsWithBonudaryPoints_ )
					{grids.push_back(FullGrid::createLinearBoundaryFullGrid(v->size(),v));}
			    else
					{grids.push_back(FullGrid::createLinearFullGrid(v->size(),v));}
				size++;
				v->pop_back();
		}
		else
		{
			for (size_t i=0;i<=sum;i++)
			{
				v->push_back(i);				
				getsums(v,dim-1,sum-i);
				v->pop_back();
			}
		}
}

void FullGridSet::getInnersums(vector<level_t> *v,size_t dim,size_t sum)
{
/** Takes recursively every possible combination of numbers which add up to sum creating a linear grid for each one(nonzero levels)*/
		if (dim==1)
		{
			    v->push_back(sum);
			    if (gridsWithBonudaryPoints_ )
					{grids.push_back(FullGrid::createLinearBoundaryFullGrid(v->size(),v));}
			    else
					{grids.push_back(FullGrid::createLinearFullGrid(v->size(),v));}
				size++;
				v->pop_back();
		}
		else
		{
			for (size_t i=1;i<=sum-1;i++)
			{
				v->push_back(i);
 				getInnersums(v,dim-1,sum-i);
				v->pop_back();
			}
		}
}

void FullGridSet::getInnersums(vector<level_t> *v,size_t dim,int sum,double* ratio)
{
/*Takes recursively every possible combination of numbers which add up to sum creating a linear grid for each one(nonzero levels)*/
		if (dim==1)
		{
				v->push_back(sum/ratio[v->size()]+1);
			    if (gridsWithBonudaryPoints_ )
					{grids.push_back(FullGrid::createLinearBoundaryFullGrid(v->size(),v));}
			    else
					{grids.push_back(FullGrid::createLinearFullGrid(v->size(),v));}
				size++;
				v->pop_back();
		}
		else{
			for (int i=0;i<=sum;i++)
			{
				v->push_back(i/ratio[v->size()]+1);
				getInnersums(v,dim-1,sum-i,ratio);
				v->pop_back();
			}
		}
}

void FullGridSet::getTrapezoidsums(vector<level_t> *v,size_t dim,size_t sum,size_t l_user)
{
	/*Takes recursively every possible combination of numbers which add up to sum creating a linear boundary grid for each one
	 * The levels of the fullgrids must be greater than l_user*/
		if (dim==1)
		{
			    v->push_back(sum);
			    if (gridsWithBonudaryPoints_ )
					{grids.push_back(FullGrid::createLinearBoundaryFullGrid(v->size(),v));}
			    else
					{grids.push_back(FullGrid::createLinearFullGrid(v->size(),v));}
				size++;
				v->pop_back();
		}
		else
		{
			for (size_t i=l_user;i<=sum-l_user;i++)
			{
				v->push_back(i);				
 				getTrapezoidsums(v,dim-1,sum-i,l_user);
				v->pop_back();
			}
		}
}

void FullGridSet::getTrapezoidsums(vector<level_t> *v,size_t dim,size_t sum,size_t* l_user)
{
	/* Takes recursively every possible combination of numbers which add up to sum creating a linear boundary grid for each one
	 * The levels of the fullgrids must be greater than l_user*/
		if (sum>=l_user[v->size()])
		{
	        if (dim==1)
			{
				    v->push_back(sum);
				    if (gridsWithBonudaryPoints_ )
						{grids.push_back(FullGrid::createLinearBoundaryFullGrid(v->size(),v));}
				    else
						{grids.push_back(FullGrid::createLinearFullGrid(v->size(),v));}
					size++;
					v->pop_back();
			}
			else
			{
				for (int i=l_user[v->size()];i<=(int)sum;i++)
				{
					v->push_back(i);
					getTrapezoidsums(v,dim-1,sum-i,l_user);
					v->pop_back();
				}
			}
		}
}

void FullGridSet::getTrapezoidsums(vector<level_t> *v,size_t dim,int sum,int l_user,double* ratio)
{
	/* Takes recursively every possible combination of numbers which add up to sum creating a linear boundary grid for each one
	 * The levels of the fullgrids must be greater than l_user*/

		if (dim==1)
		{
				v->push_back(sum/ratio[v->size()]+l_user);
			    if (gridsWithBonudaryPoints_ )
					{grids.push_back(FullGrid::createLinearBoundaryFullGrid(v->size(),v));}
			    else
					{grids.push_back(FullGrid::createLinearFullGrid(v->size(),v));}
				size++;
				v->pop_back();
		}
		else
		{
			for (int i=0;i<=sum;i++)
			{
				v->push_back(i/ratio[v->size()]+l_user);
 				getTrapezoidsums(v,dim-1,sum-i,l_user,ratio);
				v->pop_back();
			}
		}
}

void FullGridSet::getTrapezoidsums(vector<level_t> *v,size_t dim,int sum,size_t* l_user,double* ratio)
{
	/* Takes recursively every possible combination of numbers which add up to sum creating a linear boundary grid for each one
	 * The levels of the fullgrids must be greater than l_user*/
			if (dim==1)
			{
					v->push_back(sum/ratio[v->size()]+l_user[v->size()]);
				    if (gridsWithBonudaryPoints_ )
						{grids.push_back(FullGrid::createLinearBoundaryFullGrid(v->size(),v));}
				    else
						{grids.push_back(FullGrid::createLinearFullGrid(v->size(),v));}
					size++;
					v->pop_back();
			}
			else{
				for (int i=0;i<=sum;i++)
				{
					v->push_back(i/ratio[v->size()]+l_user[v->size()]);
					getTrapezoidsums(v,dim-1,sum-i,l_user,ratio);
					v->pop_back();
				}
			}
}

void FullGridSet::getSquare(vector<level_t> *v, size_t dim,size_t maxlevel)
{
    size_t small_level=maxlevel/2;
    /* Change here to the following code to take the [n/2]+1 grid as small level for odd numbers(and also change HashGenerator squareRoot method)
     * int small_level=ceil(level/2);
     * if (level%2==0) level--; */
    //we'll construct the fullgrids of level [n/2,n/2,...n,n/2,n/2] first;

	for (size_t i=0;i<dim;i++)
	{
		  v->clear();
		  for (size_t j=0;j<dim;j++)
		  {
			  if (i!=j) v->push_back(small_level);
			  else v->push_back(maxlevel);
		  }
		  // test which kind of grid to create
		  if (gridsWithBonudaryPoints_ )
		    { grids.push_back(FullGrid::createLinearBoundaryFullGrid(v->size(),v)); }
		  else
		    { grids.push_back(FullGrid::createLinearFullGrid(v->size(),v)); }
		  coefs.push_back(1);
	}
	v->clear();
	for (size_t j=0;j<dim;j++)
	{
	  v->push_back(small_level);
	}
	// test which kind of grid to create
	if (gridsWithBonudaryPoints_ )
	    { grids.push_back(FullGrid::createLinearBoundaryFullGrid(v->size(),v)); }
	else
	    { grids.push_back(FullGrid::createLinearFullGrid(v->size(),v)); }
	coefs.push_back(1.0-dim);
    size=dim+1;
}

void FullGridSet::initialize(GridStorage *storage,DataVector alpha)
{
	GridIndex* gp;
	level_t *levels=new level_t[dim];
	index_t *indexes=new index_t[dim];
	level_t *gridLevels=new level_t[dim];
	size_t i;
	int* ilev=new int[dim];
	bool good;
	int tmp,k;
	size_t j;
	double val;
	for (i=0; i < storage->size(); i++)
	{
		 gp = storage->get(i);
		 for (k=0;k<(int)dim;k++)
			  gp->get(k,levels[k],indexes[k]);
		 val=alpha[i];
		/* If the gridpoint's level is smaller than the level of the fullgrid for every dimension then the fullgrid
		 * contains the gridpoint
		 * The boolean 'good' updates for every Fullgrid as we are verifying for every coordinate the condition */
		 for (j=0;j<size;j++)
		 {
				    good=true;
				    FullGrid *fg=grids.at(j);
				    gridLevels=fg->getLevel();
				 	for (k=0;k<(int)dim;k++)
				    {
				 		ilev[k]=gridLevels[k]-levels[k];
				 		if (ilev[k]<0){good=false;break;}
				    }
				 	if (good){
				 		tmp=0;
				 		for (k=dim-1;k>=0;--k)
				 		{
				 			tmp=tmp+(tmp<<gridLevels[k])+(indexes[k]<<ilev[k]);
			 			}
				 		(*(fg))[tmp]=val;
				 	}
		   }

		/** The corresponding coordinate of a gridpoint in dimension d for level lev and index ind is
		  * new_index=ind*2^(level-lev);
		  * The position of the Fullgrid in the grids vector will be k=i1+i1*2^l1+i2*2^(l1+l2)+... */
	}
	delete[] ilev;
	delete[] levels;
	delete[] indexes;
}

void FullGridSet::deCompose(GridStorage *storage,DataVector alpha)
{

	index_t *ind=new index_t[dim];
	int *putere=new int[dim];
	GridIndex *hgi=new GridIndex(dim);
	level_t *lev;
	size_t m;
	size_t j,k;
	//div_t result;
	int aux;
	/*Takes every fullgrid from the storage and assigns a value to the corresponding element of the data vector in the fullgrid   */
	size_t startindex=grids.at(0)->startindex();
	for (size_t i=0;i<grids.size();i++)
	{
	      m=grids.at(i)->getSize();
	      lev=grids.at(i)->getLevel();
	      for (j=0;j<dim;j++) 
	  			putere[j]=grids.at(i)->length(j);
	      for (j=0;j<m;j++)		  
	      {		  
		      aux=j;	                      
		      for (k=0;k<dim;k++)
		      {
//		    	  result=div(aux,putere[k]);
//			      ind[k]=result.rem+startindex;
//		          aux=result.quot;
		    	  ind[k]=aux%putere[k]+startindex;
		    	  aux=aux/putere[k];
		      }		
   /*The level and index of the element in the hashgridstorage are computed dividing by two the index and level in the fullgrid
    *until we obtain an impair number for the index, thus obtaining the level and index in the hierarchical basis */
		      for (k=0;k<dim;k++)
		      {
			      aux=lev[k];
			      if (ind[k]!=0)
					  while (ind[k] %2==0){
						  ind[k]/=2;
						  aux--;
					  }
			      else aux=0;
		 		  hgi->push(k,aux,ind[k]);
		      }
    /*We want to obtain the position corresponding to the obtained Hashgridindex in the datavector alpha, and get the value from there*/
		      hgi->rehash();
		      aux=(*storage)[hgi];
		      (*(grids.at(i)))[j]=alpha[aux];
           }
	}
	delete[] ind;
	delete[] putere;
	delete hgi;
}

void FullGridSet::reCompose(GridStorage *storage,DataVector *alpha)
{

	GridIndex *hgi=new GridIndex(dim);

	index_t *ind=new index_t[dim];
	int *putere=new int[dim];
	level_t *lev;
	size_t m;
	size_t i,j,k;
	//div_t result;
	int aux;
	for (i=0;i<storage->size();i++)
	  (*alpha)[i]=0;
	i=0;
	size_t startindex=grids.at(0)->startindex();
	/*Takes every fullgrid from the storage and assigns a value to the corresponding element of the data vector in the fullgrid
	*We use the fact that the fullgrids are ordered by first order norm, so it is not necessary to compute the coefficients for every
	single one, just for every group with the same norm, we use the notation q=n+dim-1-norm1*/

	while (i<size)
	{
			  m=grids.at(i)->getSize();
			  lev=grids.at(i)->getLevel();
			  for (j=0;j<dim;j++)
					putere[j]=grids.at(i)->length(j);
			  for (j=0;j<m;j++)
			  {
				  aux=j;
				  for (k=0;k<dim;k++)
				  {
					//  result=div(aux,putere[k]);
					 // ind[k]=result.rem+startindex;
					 //aux=result.quot;
					  ind[k]=aux%putere[k]+startindex;
					  aux=aux/putere[k];
				  }
	   /*The level and index of the element in the hashgridstorage are computed dividing by two the index and level in the fullgrid
		*until we obtain an impair number for the index, thus obtaining the level and index in the hierarchical basis*/
				  for (k=0;k<dim;k++)
				  {
					  aux=lev[k];
					  if (ind[k]==0) aux=0;
					  else
					  while (ind[k] %2==0){
						  ind[k]/=2;
						  aux--;
					  }
					  hgi->push(k,aux,ind[k]);
				  }
		/*We want to obtain the position corresponding to the obtained Hashgridindex in the datavector alpha, and update it's value*/
				  hgi->rehash();
				  aux=(*storage)[hgi];
			      (*alpha)[aux]+=(*grids.at(i))[j]*coefs[i];
			   }
			   i++;
	}

	delete[] ind;
	delete[] putere;
	delete hgi;
}

double FullGridSet::combinedResult()
{
     double val=0;
     size_t j=0;
     /*The values obtained on fullgrids are combined after the standard combinational formula */
	 for (j=0;j<size;j++)
	    val+=grids.at(j)->val()*coefs[j];
	 return val;
}

double FullGridSet::eval(DataVector& p)
{
	double val=0.0;
	size_t j=0;
	/* Every fullgrid is evaluated and then the results are combined*/
	//std::cout << " FullGridSet::eval p[0]:" << p[0] << " , p[1]:" << p[1] << std::endl;
	for (j=0;j<size;j++)
	{
		val += grids.at(j)->eval(p) * coefs[j] ;
		//std::cout << " FullGridSet::eval j:" << j << std::endl;
		//val+=grids.at(j)->val()*coefs[j];
	}
	return val;
}

}
}

