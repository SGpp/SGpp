#include "FullGrid.hpp"
#include "LinearFullGrid.hpp"
using namespace std;
namespace sg{

    FullGrid* FullGrid::createLinearFullGrid(size_t dim, vector<level_t> *inlevel)
    {
        return new LinearFullGrid(dim,inlevel);
    }

    FullGrid* FullGrid::createLinearBoundaryFullGrid(size_t dim, vector<level_t> *inlevel)
    {
        return new FullGrid(dim,inlevel);
    }

   void FullGrid::getCoords(size_t index,DataVector &v)
    {
        int ind=0;
        int aux=index;
        /**
         *The number of gridpoints on every direction is Ni:=2^level[i]+1, so the index of a gridpoint in a given direction
         *can be computed by the formula ind(i)=[index/(N0*N1...        *Ni-1)] mod Ni, where index is the position of the gridpoint
         *in the vector.
         *The coordinate of a gridpoint in direction 'i' is ind(i)*2^(-level[i]), considering the grid of dimension 1.0
         */
        for (size_t i=0;i<dim;i++){
               ind=aux%((int)powOfTwo[level[i]]+1);
	           aux=aux/((int)powOfTwo[level[i]]+1);
               if (ind==0) v[i]=0;
               else
               v[i]=((double)ind)*negPowOfTwo[level[i]];
        }
        if (FullGrid::boundingBox!=0)
        {
        	 for (size_t i=0;i<dim;i++)
        		 v[i]=v[i]*FullGrid::boundingBox->getIntervalWidth(i)+FullGrid::boundingBox->getIntervalOffset(i);
        }
    }

   std::string FullGrid::getCoordsString(size_t index)
    {
    	std::stringstream return_stream;

    	// switch on scientific notation:
		//return_stream << std::scientific;
        DataVector coords(dim);
        getCoords(index,coords);
        //return_stream.precision(15);
    	for(size_t i = 0; i < dim; i++)
    	{
    		return_stream << std:: scientific<<coords[i];    		
    		if (i < dim-1)
    		{
    			return_stream << " ";
    		}
    	}

    	return return_stream.str();
    }
}
