#ifndef FDMSTENCILS_HPP_
#define FDMSTENCILS_HPP_



template<typename GIT, typename STORAGE, typename DATA>
struct LABSHierarchSurplusStencil
{
	void operator() (STORAGE& storage, GIT &index, int d, DATA &interm, DATA &alpha)
	{
        	index.set(d, 0, 0);
        	double lb = interm[storage[&index]];
			alpha[storage[&index]] = lb;

        	index.set(d, 0, 1);
        	double rb = interm[storage[&index]];
        	alpha[storage[&index]] = rb;
			
			index.set(d, 1, 1);
			hierarchSurplus(storage, lb, rb, index, d, interm, alpha);		
	}
	
	void hierarchSurplus(STORAGE& storage, double lb, double rb, GIT &index, int d, DATA &interm, DATA &alpha)
	{
		double interpolated = (lb + rb)/2.0;
		
		typename STORAGE::grid_iterator iter = storage.find(&index);
		double value = interm[iter->second];
		alpha[iter->second] = value - interpolated;
		
		
		typename GIT::index_type srcIndex;
        typename GIT::level_type srcLevel;
        
        index.get(d, srcLevel, srcIndex);
        
        index.set(d, srcLevel+1, srcIndex*2 - 1);
        if(storage.has_key(&index))
        {
        	hierarchSurplus(storage, lb, value, index, d, interm, alpha);
        }

        index.set(d, srcLevel+1, srcIndex*2 + 1);
        if(storage.has_key(&index))
        {
        	hierarchSurplus(storage, value, rb, index, d, interm, alpha);
        }

        index.push(d, srcLevel, srcIndex);		
	}
		
};

template<typename GIT, typename STORAGE, typename DATA>
struct LABSDehierarchSurplusStencil
{
	void operator() (STORAGE& storage, GIT &index, int d, DATA &interm, DATA &alpha)
	{
        	index.set(d, 0, 0);
        	double lb = interm[storage[&index]];
			alpha[storage[&index]] = lb;

        	index.set(d, 0, 1);
        	double rb = interm[storage[&index]];
        	alpha[storage[&index]] = rb;
			
			index.set(d, 1, 1);
			hierarchSurplus(storage, lb, rb, index, d, interm, alpha);		
	}
	
	void hierarchSurplus(STORAGE& storage, double lb, double rb, GIT &index, int d, DATA &interm, DATA &alpha)
	{
		double interpolated = (lb + rb)/2.0;
		
		typename STORAGE::grid_iterator iter = storage.find(&index);
		double value = interm[iter->second] + interpolated;
		alpha[iter->second] = value;
		
		
		typename GIT::index_type srcIndex;
        typename GIT::level_type srcLevel;
        
        index.get(d, srcLevel, srcIndex);
        
        index.set(d, srcLevel+1, srcIndex*2 - 1);
        if(storage.has_key(&index))
        {
        	hierarchSurplus(storage, lb, value, index, d, interm, alpha);
        }

        index.set(d, srcLevel+1, srcIndex*2 + 1);
        if(storage.has_key(&index))
        {
        	hierarchSurplus(storage, value, rb, index, d, interm, alpha);
        }

        index.push(d, srcLevel, srcIndex);		
	}
		
};


#endif /*FDMSTENCILS_HPP_*/
