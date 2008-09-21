#ifndef LABSSURPLUS_HPP_
#define LABSSURPLUS_HPP_


template<typename GIT, typename STORAGE, typename DATA>
class LABSSurplus
{
public:
	LABSSurplus(STORAGE &storage) : storage(storage)
	{
	}

	void hierarchSurplus(DATA &alpha)
	{
		for(int i = 0; i < GIT::dim; i++)
		{
			hierarchSurplus1D(i, alpha, alpha);
		}	
	}

	void hierarchSurplus1D(int dim, DATA &interm, DATA &alpha)
	{
		GIT index;
		
        for(int d = 0; d < GIT::dim; d++)
        {
            index.push(d,0,0);
        }
        index.rehash();
        
        if(GIT::dim == 1)
        {
        	index.set(0, 0, 0);
        	double lb = interm[storage[&index]];
			alpha[storage[&index]] = lb;

        	index.set(0, 0, 1);
        	double rb = interm[storage[&index]];
        	alpha[storage[&index]] = rb;
        	
        	index.set(0, 1, 1);
        	hierarchSurplus1D(lb, rb, index, 0, interm, alpha);
        }
        else
        {
			if(dim != 0)
			{
				hierarchSurplus(0, index, dim, interm, alpha);
			} 
			else
			{
				hierarchSurplus(1, index, dim, interm, alpha);
			}
        }
	}

protected:

	void hierarchSurplus(int recD, GIT &index, int d, DATA &interm, DATA &alpha)
	{
		typename GIT::index_type srcIndex;
        typename GIT::level_type srcLevel;

		int nextD = recD + 1;
		if(nextD == d)
		{
			nextD++;
		}

		index.get(recD, srcLevel, srcIndex);

		bool resetToZero = false;

		if(srcLevel == 0)
		{
			index.set(recD, 0, 0);
			hierarchSurplusRec(nextD, index, d, interm, alpha);

			index.set(recD, 0, 1);
			hierarchSurplusRec(nextD, index, d, interm, alpha);
			
			srcLevel = 1;
			srcIndex = 1;
			resetToZero = true;
		}

		
		
		index.set(recD, srcLevel+1, srcIndex*2 - 1);
		if(storage.has_key(&index))
		{
			hierarchSurplus(recD, index, d, interm, alpha);
		}
		
		index.set(recD, srcLevel+1, srcIndex*2 + 1);
		if(storage.has_key(&index))
		{
			hierarchSurplus(recD, index, d, interm, alpha);
		}
		
		// TODO
		//index.set(recD, 0, 0);
		hierarchSurplusRec(nextD, index, d, interm, alpha);		
	}

	void hierarchSurplusRec(int recD, GIT &index, int d, DATA &interm, DATA &alpha)
	{
		if(recD >= GIT::dim)
		{
        	index.set(d, 0, 0);
        	double lb = interm[storage[&index]];
			alpha[storage[&index]] = lb;

        	index.set(d, 0, 1);
        	double rb = interm[storage[&index]];
        	alpha[storage[&index]] = rb;
			
			index.set(d, 1, 1);
			hierarchSurplus1D(lb, rb, index, d, interm, alpha);
		}
		else
		{
			hierarchSurplus(recD, index, d, interm, alpha);
		}
	}

	void hierarchSurplus1D(double lb, double rb, GIT &index, int d, DATA &interm, DATA &alpha)
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
        	hierarchSurplus1D(lb, value, index, d, interm, alpha);
        }

        index.set(d, srcLevel+1, srcIndex*2 + 1);
        if(storage.has_key(&index))
        {
        	hierarchSurplus1D(value, rb, index, d, interm, alpha);
        }

        index.push(d, srcLevel, srcIndex);		
	}


private:
	STORAGE &storage;
		
};


#endif /*LABSSURPLUS_HPP_*/
